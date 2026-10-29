//! Filesystem publication policy, separate from numerical command execution.
use std::fs::{self, File, OpenOptions};
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::{SystemTime, UNIX_EPOCH};

use serde::Serialize;
use sha2::{Digest, Sha256};

use crate::errors::DeseqError;

static TEMP_ID: AtomicU64 = AtomicU64::new(0);

struct TemporaryPath(PathBuf);
impl Drop for TemporaryPath {
    fn drop(&mut self) {
        if self.0.is_dir() {
            let _ = fs::remove_dir_all(&self.0);
        } else {
            let _ = fs::remove_file(&self.0);
        }
    }
}

struct StagedFile {
    temporary: TemporaryPath,
    destination: PathBuf,
}

pub(super) struct OutputTransaction {
    files: Vec<StagedFile>,
    marker: StagedFile,
}

#[derive(Serialize)]
struct CompletedFile {
    path: PathBuf,
    bytes: u64,
    sha256: String,
}

#[derive(Serialize)]
struct CompletionManifest {
    schema_version: u32,
    generation: String,
    rsdeseq2_version: &'static str,
    files: Vec<CompletedFile>,
}

fn invalid(reason: impl Into<String>) -> DeseqError {
    DeseqError::InvalidOptions {
        reason: reason.into(),
    }
}

fn destination(path: &Path) -> Result<PathBuf, DeseqError> {
    if path.exists() {
        let resolved = path.canonicalize()?;
        if !resolved.is_file() {
            return Err(invalid(format!(
                "output destination is not a file: {}",
                path.display()
            )));
        }
        Ok(resolved)
    } else {
        // Canonicalizing the existing parent also resolves relative paths and
        // directory symlinks. Missing parents fail before any fitting or writes.
        let parent = path
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
            .unwrap_or(Path::new("."));
        let name = path
            .file_name()
            .ok_or_else(|| invalid("output requires a file name"))?;
        Ok(parent.canonicalize()?.join(name))
    }
}

fn unique_temporary(parent: &Path, directory: bool) -> Result<TemporaryPath, DeseqError> {
    for _ in 0..128 {
        let id = TEMP_ID.fetch_add(1, Ordering::Relaxed);
        let path = parent.join(format!(".rsdeseq2-{}-{id}.tmp", std::process::id()));
        let created = if directory {
            fs::create_dir(&path)
        } else {
            OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(&path)
                .map(drop)
        };
        match created {
            Ok(()) => return Ok(TemporaryPath(path)),
            Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => continue,
            Err(error) => return Err(error.into()),
        }
    }
    Err(invalid("unable to reserve an output staging path"))
}

fn aliases(a: &Path, b: &Path) -> Result<bool, DeseqError> {
    Ok(a == b || (a.exists() && b.exists() && same_file::is_same_file(a, b)?))
}

fn hash_file(path: &Path) -> Result<(u64, String), DeseqError> {
    let mut file = File::open(path)?;
    let mut hash = Sha256::new();
    let mut bytes = 0;
    let mut buffer = [0_u8; 64 * 1024];
    loop {
        let count = file.read(&mut buffer)?;
        if count == 0 {
            break;
        }
        hash.update(&buffer[..count]);
        bytes += count as u64;
    }
    let digest = hash.finalize();
    let mut encoded = String::with_capacity(digest.len() * 2);
    const HEX: &[u8; 16] = b"0123456789abcdef";
    for byte in digest {
        encoded.push(HEX[(byte >> 4) as usize] as char);
        encoded.push(HEX[(byte & 0x0f) as usize] as char);
    }
    Ok((bytes, encoded))
}

impl OutputTransaction {
    pub(super) fn prepare(
        inputs: &[&mut PathBuf],
        outputs: &mut [&mut PathBuf],
    ) -> Result<Self, DeseqError> {
        let inputs = inputs
            .iter()
            .map(|path| path.canonicalize())
            .collect::<Result<Vec<_>, _>>()?;
        let mut destinations = outputs
            .iter()
            .map(|path| destination(path))
            .collect::<Result<Vec<_>, _>>()?;
        // The primary output is last in the explicit path inventory.
        let main = destinations
            .last()
            .ok_or_else(|| invalid("command has no output path"))?;
        let mut marker_name = main.file_name().unwrap().to_os_string();
        marker_name.push(".complete.json");
        destinations.push(destination(&main.with_file_name(marker_name))?);
        let mut probes: Vec<(PathBuf, TemporaryPath)> = Vec::new();
        for (index, output) in destinations.iter().enumerate() {
            for input in &inputs {
                if aliases(output, input)? {
                    return Err(invalid(format!(
                        "input/output path collision: {}",
                        output.display()
                    )));
                }
            }
            for previous in &destinations[..index] {
                if aliases(output, previous)? {
                    return Err(invalid(format!(
                        "output path collision: {}",
                        output.display()
                    )));
                }
            }
            let parent = output.parent().unwrap();
            let mut probe_index = None;
            for (i, (existing_parent, _)) in probes.iter().enumerate() {
                if aliases(parent, existing_parent)? {
                    probe_index = Some(i);
                    break;
                }
            }
            let i = if let Some(i) = probe_index {
                i
            } else {
                probes.push((parent.to_owned(), unique_temporary(parent, true)?));
                probes.len() - 1
            };
            // Ask the destination filesystem itself about case/Unicode aliases;
            // no assumptions about the host's filesystem or string case folding.
            let probe = probes[i].1.0.join(output.file_name().unwrap());
            match OpenOptions::new().write(true).create_new(true).open(probe) {
                Ok(file) => drop(file),
                Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => {
                    return Err(invalid(format!(
                        "output path collision on destination filesystem: {}",
                        output.display()
                    )));
                }
                Err(error) => return Err(error.into()),
            }
        }
        let mut staged = Vec::new();
        for destination in destinations {
            staged.push(StagedFile {
                temporary: unique_temporary(destination.parent().unwrap(), false)?,
                destination,
            });
        }
        let marker = staged.pop().unwrap();
        for (output, file) in outputs.iter_mut().zip(&staged) {
            **output = file.temporary.0.clone();
        }
        Ok(Self {
            files: staged,
            marker,
        })
    }

    pub(super) fn publish(self) -> Result<(), DeseqError> {
        self.publish_with(|_, temporary, destination| fs::rename(temporary, destination))
    }

    fn publish_with(
        self,
        mut rename: impl FnMut(usize, &Path, &Path) -> std::io::Result<()>,
    ) -> Result<(), DeseqError> {
        let mut completed = Vec::new();
        for file in &self.files {
            OpenOptions::new()
                .write(true)
                .open(&file.temporary.0)?
                .sync_all()?;
            let (bytes, sha256) = hash_file(&file.temporary.0)?;
            if bytes == 0 {
                return Err(invalid("requested output was not written"));
            }
            completed.push(CompletedFile {
                path: file.destination.clone(),
                bytes,
                sha256,
            });
        }
        let generation = format!(
            "{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .map_err(|e| invalid(e.to_string()))?
                .as_nanos()
        );
        let manifest = CompletionManifest {
            schema_version: 1,
            generation,
            rsdeseq2_version: env!("CARGO_PKG_VERSION"),
            files: completed,
        };
        let mut marker = File::create(&self.marker.temporary.0)?;
        serde_json::to_writer_pretty(&mut marker, &manifest).map_err(|e| invalid(e.to_string()))?;
        marker.write_all(b"\n")?;
        marker.sync_all()?;
        drop(marker);
        // Invalidate the previous generation only after all staged outputs have
        // completed. Separate renames are not a multi-file atomic transaction.
        match fs::remove_file(&self.marker.destination) {
            Ok(()) => {}
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => {}
            Err(error) => return Err(error.into()),
        }
        for (index, file) in self.files.iter().enumerate() {
            rename(index, &file.temporary.0, &file.destination)?;
        }
        rename(
            self.files.len(),
            &self.marker.temporary.0,
            &self.marker.destination,
        )?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn publication_failure_invalidates_the_marker_and_cleans_temporary_files() {
        let root = unique_temporary(&std::env::temp_dir(), true).unwrap();
        let first_final = root.0.join("first.tsv");
        let second_final = root.0.join("second.tsv");
        let marker = root.0.join("second.tsv.complete.json");
        fs::write(&first_final, "old first").unwrap();
        fs::write(&second_final, "old second").unwrap();
        fs::write(&marker, "previous generation").unwrap();
        let mut first = first_final.clone();
        let mut second = second_final.clone();
        let transaction = OutputTransaction::prepare(&[], &mut [&mut first, &mut second]).unwrap();
        fs::write(&first, "new first").unwrap();
        fs::write(&second, "new second").unwrap();
        let error = transaction
            .publish_with(|index, from, to| {
                if index == 1 {
                    Err(std::io::Error::other("injected rename failure"))
                } else {
                    fs::rename(from, to)
                }
            })
            .unwrap_err();
        assert!(error.to_string().contains("injected rename failure"));
        assert_eq!(fs::read_to_string(first_final).unwrap(), "new first");
        assert_eq!(fs::read_to_string(second_final).unwrap(), "old second");
        assert!(!marker.exists());
        assert!(!first.exists() && !second.exists());
        assert_eq!(fs::read_dir(&root.0).unwrap().count(), 2);
    }

    #[test]
    fn complete_manifest_hashes_all_successful_outputs() {
        let root = unique_temporary(&std::env::temp_dir(), true).unwrap();
        let final_path = root.0.join("result.tsv");
        let mut output = final_path.clone();
        let transaction = OutputTransaction::prepare(&[], &mut [&mut output]).unwrap();
        fs::write(&output, "abc").unwrap();
        transaction.publish().unwrap();
        let manifest: serde_json::Value =
            serde_json::from_reader(File::open(root.0.join("result.tsv.complete.json")).unwrap())
                .unwrap();
        assert_eq!(manifest["schema_version"], 1);
        assert_eq!(manifest["files"][0]["bytes"], 3);
        assert_eq!(
            manifest["files"][0]["path"],
            final_path.canonicalize().unwrap().to_str().unwrap()
        );
        assert_eq!(
            manifest["files"][0]["sha256"],
            "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
        );
        assert_eq!(fs::read_to_string(final_path).unwrap(), "abc");
        assert!(!output.exists());
        assert_eq!(fs::read_dir(&root.0).unwrap().count(), 2);
    }

    #[test]
    fn unwritten_output_preserves_the_previous_generation() {
        let root = unique_temporary(&std::env::temp_dir(), true).unwrap();
        let final_path = root.0.join("result.tsv");
        let marker = root.0.join("result.tsv.complete.json");
        fs::write(&final_path, "old result").unwrap();
        fs::write(&marker, "old marker").unwrap();
        let mut output = final_path.clone();
        let transaction = OutputTransaction::prepare(&[], &mut [&mut output]).unwrap();
        assert!(
            transaction
                .publish()
                .unwrap_err()
                .to_string()
                .contains("not written")
        );
        assert_eq!(fs::read_to_string(final_path).unwrap(), "old result");
        assert_eq!(fs::read_to_string(marker).unwrap(), "old marker");
        assert!(!output.exists());
    }

    #[test]
    fn preflight_detects_existing_identity_and_filesystem_name_aliases() {
        let root = unique_temporary(&std::env::temp_dir(), true).unwrap();
        let mut input = root.0.join("counts.tsv");
        fs::write(&input, "counts").unwrap();
        let mut hardlink = root.0.join("hardlink.tsv");
        fs::hard_link(&input, &hardlink).unwrap();
        assert!(OutputTransaction::prepare(&[&mut input], &mut [&mut hardlink]).is_err());
        #[cfg(unix)]
        {
            let mut symlink = root.0.join("symlink.tsv");
            std::os::unix::fs::symlink(&input, &symlink).unwrap();
            assert!(OutputTransaction::prepare(&[&mut input], &mut [&mut symlink]).is_err());
        }
        let mut first = root.0.join("Result.tsv");
        let mut second = root.0.join("result.tsv");
        let case_probe = root.0.join("CaseProbe");
        fs::write(case_probe, "probe").unwrap();
        let insensitive = root.0.join("caseprobe").exists();
        let result = OutputTransaction::prepare(&[], &mut [&mut first, &mut second]);
        assert_eq!(result.is_err(), insensitive);
        assert_eq!(fs::read_to_string(input).unwrap(), "counts");
    }
}
