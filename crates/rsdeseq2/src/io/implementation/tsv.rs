use crate::errors::DeseqError;
use csv_core::{ReadRecordResult, WriteResult};
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::path::Path;

#[derive(Default)]
pub(super) struct Record<'a> {
    text: &'a str,
    ends: &'a [usize],
}

impl<'a> Record<'a> {
    pub(super) fn len(&self) -> usize {
        self.ends.len()
    }

    pub(super) fn get(&self, index: usize) -> Option<&'a str> {
        let &end = self.ends.get(index)?;
        let start = if index == 0 { 0 } else { self.ends[index - 1] };
        Some(&self.text[start..end])
    }

    pub(super) fn iter(&self) -> impl Iterator<Item = &'a str> + '_ {
        std::iter::once(0)
            .chain(self.ends.iter().copied())
            .zip(self.ends.iter().copied())
            .map(|(start, end)| &self.text[start..end])
    }
}

pub(super) struct TsvReader<R> {
    source: BufReader<R>,
    parser: csv_core::Reader,
    data: Vec<u8>,
    fields: Vec<usize>,
    width: Option<usize>,
}

impl<R: Read> TsvReader<R> {
    pub(super) fn new(source: R) -> Self {
        Self {
            source: BufReader::new(source),
            parser: csv_core::ReaderBuilder::new().delimiter(b'\t').build(),
            data: vec![0; 512],
            fields: vec![0; 32],
            width: None,
        }
    }

    pub(super) fn next_record(&mut self) -> Result<Option<Record<'_>>, DeseqError> {
        let line = self.parser.line();
        let mut nbytes = 0;
        let mut nfields = 0;
        loop {
            let (status, consumed, produced, fields) = self.parser.read_record(
                self.source.fill_buf()?,
                &mut self.data[nbytes..],
                &mut self.fields[nfields..],
            );
            self.source.consume(consumed);
            nbytes += produced;
            nfields += fields;
            match status {
                ReadRecordResult::InputEmpty => {}
                ReadRecordResult::OutputFull => self.data.resize(self.data.len() * 2, 0),
                ReadRecordResult::OutputEndsFull => self.fields.resize(self.fields.len() * 2, 0),
                ReadRecordResult::Record => break,
                ReadRecordResult::End => return Ok(None),
            }
        }

        let expected = *self.width.get_or_insert(nfields);
        if expected != nfields {
            return Err(DeseqError::InvalidTableWidth {
                line,
                expected,
                actual: nfields,
            });
        }
        let text = std::str::from_utf8(&self.data[..nbytes])
            .map_err(|_| DeseqError::InvalidTableEncoding { line })?;
        let ends = &self.fields[..nfields];
        if ends.iter().any(|&end| !text.is_char_boundary(end)) {
            return Err(DeseqError::InvalidTableEncoding { line });
        }
        Ok(Some(Record { text, ends }))
    }
}

pub(crate) struct TsvWriter {
    destination: BufWriter<File>,
    writer: csv_core::Writer,
    buffer: [u8; 8192],
    width: Option<usize>,
}

impl TsvWriter {
    pub(crate) fn from_path(path: impl AsRef<Path>) -> io::Result<Self> {
        Ok(Self {
            destination: BufWriter::new(File::create(path)?),
            writer: csv_core::WriterBuilder::new().delimiter(b'\t').build(),
            buffer: [0; 8192],
            width: None,
        })
    }

    pub(crate) fn write_record(
        &mut self,
        record: impl IntoIterator<Item = impl AsRef<[u8]>>,
    ) -> io::Result<()> {
        let mut width = 0;
        for field in record {
            if width != 0 {
                let (_, written) = self.writer.delimiter(&mut self.buffer);
                self.destination.write_all(&self.buffer[..written])?;
            }
            let mut input = field.as_ref();
            loop {
                let (status, consumed, written) = self.writer.field(input, &mut self.buffer);
                self.destination.write_all(&self.buffer[..written])?;
                if status == WriteResult::InputEmpty {
                    break;
                }
                input = &input[consumed..];
            }
            width += 1;
        }
        if width != *self.width.get_or_insert(width) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "table output record width differs from its header",
            ));
        }
        let (_, written) = self.writer.terminator(&mut self.buffer);
        self.destination.write_all(&self.buffer[..written])
    }

    pub(crate) fn flush(&mut self) -> io::Result<()> {
        self.destination.flush()
    }
}
