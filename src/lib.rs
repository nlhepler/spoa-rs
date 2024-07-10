use spoa_sys::ffi;
use std::ffi::{c_char, CStr, CString};

pub use ffi::AlignmentType;

/// An opaque struct containing an alignment of a sequence to a `Graph`
pub struct Alignment(cxx::UniquePtr<ffi::Alignment>);

/// An opaque struct for aligning sequences to a `Graph`
pub struct AlignmentEngine(cxx::UniquePtr<ffi::AlignmentEngine>);

impl AlignmentEngine {
    /// Construct a new `AlignmentEngine`, given an `AlignmentType`, a match score (`m`), a
    /// mismatch score (`n`), a gap open score (`g`), a gap extension score (`e`), a second gap
    /// open score (`q`), and a second gap extension score (`c`).
    pub fn new(typ: AlignmentType, m: i8, n: i8, g: i8, e: i8, q: i8, c: i8) -> Self {
        AlignmentEngine(ffi::create_alignment_engine(typ, m, n, g, e, q, c))
    }

    /// Align a `sequence` to a `graph`, returning an `Alignment`
    pub fn align(&mut self, sequence: &CStr, graph: &Graph) -> Alignment {
        let sequence_len = u32::try_from(sequence.to_bytes().len()).unwrap();
        Alignment(unsafe {
            ffi::align(
                self.0.pin_mut(),
                sequence.as_ptr(),
                sequence_len,
                graph.0.as_ref().unwrap(),
            )
        })
    }

    /// Align a `sequence` to a `graph`, returning an `Alignment`, from a byte slice
    /// of sequence and quality lengths.
    pub fn align_from_bytes(&mut self, sequence: &[u8], graph: &Graph) -> Alignment {
        let sequence_len: u32 = sequence
            .len()
            .try_into()
            .expect("Sequence length does not fit into a u32");
        Alignment(unsafe {
            ffi::align(
                self.0.pin_mut(),
                sequence.as_ptr() as *const c_char,
                sequence_len,
                graph.0.as_ref().unwrap(),
            )
        })
    }
}

/// An opaque struct for the partial order alignment graph
pub struct Graph(cxx::UniquePtr<ffi::Graph>);

impl Graph {
    /// Construct a new, empty `Graph`.
    pub fn new() -> Self {
        Graph(ffi::create_graph())
    }

    /// Add an `Alignment` to a `Graph`, along with its `sequence` and per-base `quality` scores.
    ///
    /// The `alignment` should be derived from calling `AlignmentEngine::align` with the
    /// `sequence`.
    pub fn add_alignment(&mut self, alignment: &Alignment, sequence: &CStr, quality: &CStr) {
        let sequence_len = u32::try_from(sequence.to_bytes().len()).unwrap();
        let quality_len = u32::try_from(quality.to_bytes().len()).unwrap();
        assert!(sequence_len == quality_len);
        unsafe {
            ffi::add_alignment(
                self.0.pin_mut(),
                alignment.0.as_ref().unwrap(),
                sequence.as_ptr(),
                sequence_len,
                quality.as_ptr(),
                quality_len,
            )
        }
    }

    /// Add an `Alignment` to a `Graph`, along with its `sequence` and per-base `quality` scores, from a byte slice of
    /// sequence and quality lengths.
    ///
    /// The `alignment` should be derived from calling `AlignmentEngine::align` with the
    /// `sequence`.
    pub fn add_alignment_from_bytes(
        &mut self,
        alignment: &Alignment,
        sequence: &[u8],
        quality: &[u8],
    ) {
        let sequence_len = sequence
            .len()
            .try_into()
            .expect("Sequence length does not fit into a u32");
        let quality_len = quality
            .len()
            .try_into()
            .expect("Quality length does not fit into a u32");

        assert!(sequence_len == quality_len);

        unsafe {
            ffi::add_alignment(
                self.0.pin_mut(),
                alignment.0.as_ref().unwrap(),
                sequence.as_ptr() as *const c_char,
                sequence_len,
                quality.as_ptr() as *const c_char,
                quality_len,
            )
        }
    }

    /// Generate a consenus sequence from the partial order `Graph`.
    pub fn consensus(&mut self) -> CString {
        let mut buf = Vec::from(
            ffi::generate_consensus(self.0.pin_mut())
                .as_ref()
                .unwrap()
                .as_bytes(),
        );
        buf.push(0);
        CString::from_vec_with_nul(buf).unwrap()
    }

    /// Generate a multiple sequence alignment for all sequences added to the `Graph`.
    ///
    /// If `include_consensus` is provided, the consensus sequence is provided, also aligned, at
    /// the end.
    pub fn multiple_sequence_alignment(&mut self, include_consensus: bool) -> Vec<CString> {
        let msa = ffi::generate_multiple_sequence_alignment(self.0.pin_mut(), include_consensus);
        let mut result = vec![];
        for aln in msa.as_ref().unwrap().iter() {
            let mut buf = Vec::from(aln.as_bytes());
            buf.push(0);
            result.push(CString::from_vec_with_nul(buf).unwrap());
        }
        result
    }
}

impl Default for Graph {
    fn default() -> Self {
        Graph::new()
    }
}

#[cfg(test)]
mod fastq;

#[cfg(test)]
mod tests {
    use super::*;

    // thanks to pjedge for this small
    const SMALL_SEQS: &[&CStr] = unsafe {
        &[
            CStr::from_bytes_with_nul_unchecked(b"ATTGCCCGTT\0"),
            CStr::from_bytes_with_nul_unchecked(b"AATGCCGTT\0"),
            CStr::from_bytes_with_nul_unchecked(b"AATGCCCGAT\0"),
            CStr::from_bytes_with_nul_unchecked(b"AACGCCCGTC\0"),
            CStr::from_bytes_with_nul_unchecked(b"AGTGCTCGTT\0"),
            CStr::from_bytes_with_nul_unchecked(b"AATGCTCGTT\0"),
        ]
    };

    // same sequences as SMALL_SEQS, just as a &[u8] byte slice
    const SMALL_SEQS_BYTES: &[&[u8]] = &[
        "ATTGCCCGTT".as_bytes(),
        "AATGCCGTT".as_bytes(),
        "AATGCCCGAT".as_bytes(),
        "AACGCCCGTC".as_bytes(),
        "AGTGCTCGTT".as_bytes(),
        "AATGCTCGTT".as_bytes(),
    ];

    #[test]
    fn consensus_small() {
        let mut eng = AlignmentEngine::new(AlignmentType::kNW, 5, -4, -3, -1, -3, -1);
        let mut graph = Graph::new();

        for seq in SMALL_SEQS {
            let qual = {
                let mut qual = vec![34u8; seq.to_bytes().len()];
                qual.push(0);
                CString::from_vec_with_nul(qual).unwrap()
            };
            let aln = eng.align(seq, &graph);
            graph.add_alignment(&aln, seq, &qual);
        }

        let consensus = graph.consensus();
        assert!(consensus.as_bytes() == b"AATGCCCGTT");
    }

    #[test]
    fn consensus_small_bytes() {
        let mut eng = AlignmentEngine::new(AlignmentType::kNW, 5, -4, -3, -1, -3, -1);
        let mut graph = Graph::new();

        for seq in SMALL_SEQS_BYTES {
            let qual = vec![34u8; seq.len()];

            let aln = eng.align_from_bytes(seq, &graph);
            graph.add_alignment_from_bytes(&aln, seq, &qual);
        }

        let consensus = graph.consensus();
        assert!(consensus.as_bytes() == b"AATGCCCGTT");
    }

    #[test]
    fn consensus_spoa() {
        // test borrowed from spoa itself, obviously
        let file = std::fs::File::open("spoa-sys/spoa/test/data/sample.fastq.gz").unwrap();
        let gz = flate2::read::GzDecoder::new(file);
        let reader = fastq::FastqReader::new(std::io::BufReader::new(gz));

        let mut eng = AlignmentEngine::new(AlignmentType::kSW, 5, -4, -8, -6, -8, -6);
        let mut graph = Graph::new();

        for record in reader {
            let record = record.unwrap();
            let aln = eng.align(&record.seq, &graph);
            graph.add_alignment(&aln, &record.seq, &record.qual);
        }

        let consensus = graph.consensus();
        let expected = b"AATGATGCGCTTTGTTGGCGCGGTGGCTTGATGCAGGGGCTAATCGACCTCTGGCAACCACTTTTCCATGACAGGAGTTGAATATGGCATTCAGTAATCCCTTCGATGATCCGCAGGGAGCGTTTTACATATTGCGCAATGCGCAGGGGCAATTCAGTCTGTGGCCGCAACAATGCGTCTTACCGGCAGGCTGGGACATTGTGTGTCAGCCGCAGTCACAGGCGTCCTGCCAGCAGTGGCTGGAAGCCCACTGGCGTACTCTGACACCGACGAATTTTACCCAGTTGCAGGAGGCACAATGAGCCAGCATTTACCTTTGGTCGCCGCACAGCCCGGCATCTGGATGGCAGAAAAACTGTCAGAATTACCCTCCGCCTGGAGCGTGGCGCATTACGTTGAGTTAACCGGAGAGGTTGATTCGCCATTACTGGCCCGCGCGGTGGTTGCCGGACTAGCGCAAGCAGATACGCTTTACACGCGCAACCAAGGATTTCGG";
        assert!(consensus.as_bytes() == expected);
    }

    #[test]
    fn consensus_spoa_bytes() {
        // test borrowed from spoa itself, obviously
        let file = std::fs::File::open("spoa-sys/spoa/test/data/sample.fastq.gz").unwrap();
        let gz = flate2::read::GzDecoder::new(file);
        let reader = fastq::FastqReader::new(std::io::BufReader::new(gz));

        let mut eng = AlignmentEngine::new(AlignmentType::kSW, 5, -4, -8, -6, -8, -6);
        let mut graph = Graph::new();

        for record in reader {
            let record = record.unwrap();

            // fastq::FastqRecord uses CString types - we can
            // convert this to a byte array.
            let aln = eng.align_from_bytes(record.seq.as_bytes(), &graph);
            graph.add_alignment_from_bytes(&aln, record.seq.as_bytes(), record.qual.as_bytes());
        }

        let consensus = graph.consensus();
        let expected = b"AATGATGCGCTTTGTTGGCGCGGTGGCTTGATGCAGGGGCTAATCGACCTCTGGCAACCACTTTTCCATGACAGGAGTTGAATATGGCATTCAGTAATCCCTTCGATGATCCGCAGGGAGCGTTTTACATATTGCGCAATGCGCAGGGGCAATTCAGTCTGTGGCCGCAACAATGCGTCTTACCGGCAGGCTGGGACATTGTGTGTCAGCCGCAGTCACAGGCGTCCTGCCAGCAGTGGCTGGAAGCCCACTGGCGTACTCTGACACCGACGAATTTTACCCAGTTGCAGGAGGCACAATGAGCCAGCATTTACCTTTGGTCGCCGCACAGCCCGGCATCTGGATGGCAGAAAAACTGTCAGAATTACCCTCCGCCTGGAGCGTGGCGCATTACGTTGAGTTAACCGGAGAGGTTGATTCGCCATTACTGGCCCGCGCGGTGGTTGCCGGACTAGCGCAAGCAGATACGCTTTACACGCGCAACCAAGGATTTCGG";
        assert!(consensus.as_bytes() == expected);
    }

    #[test]
    fn msa_small() {
        let mut eng = AlignmentEngine::new(AlignmentType::kNW, 5, -4, -3, -1, -3, -1);
        let mut graph = Graph::new();

        for seq in SMALL_SEQS {
            let qual = {
                let mut qual = vec![34u8; seq.to_bytes().len()];
                qual.push(0);
                CString::from_vec_with_nul(qual).unwrap()
            };
            let aln = eng.align(seq, &graph);
            graph.add_alignment(&aln, seq, &qual);
        }

        let msa = graph.multiple_sequence_alignment(true);
        assert!(msa.len() == 7);
        assert!(msa[0].as_bytes() == b"ATTGCC-CGTT");
        assert!(msa[1].as_bytes() == b"AATG-C-CGTT");
        assert!(msa[2].as_bytes() == b"AATGCC-CGAT");
        assert!(msa[3].as_bytes() == b"AACGCC-CGTC");
        assert!(msa[4].as_bytes() == b"AGTG-CTCGTT");
        assert!(msa[5].as_bytes() == b"AATG-CTCGTT");
        assert!(msa[6].as_bytes() == b"AATGCC-CGTT");
    }

    #[test]
    fn msa_small_bytes() {
        let mut eng = AlignmentEngine::new(AlignmentType::kNW, 5, -4, -3, -1, -3, -1);
        let mut graph = Graph::new();

        for seq in SMALL_SEQS_BYTES {
            let qual = vec![34u8; seq.len()];

            let aln = eng.align_from_bytes(seq, &graph);
            graph.add_alignment_from_bytes(&aln, seq, &qual);
        }

        let msa = graph.multiple_sequence_alignment(true);
        assert!(msa.len() == 7);
        assert!(msa[0].as_bytes() == b"ATTGCC-CGTT");
        assert!(msa[1].as_bytes() == b"AATG-C-CGTT");
        assert!(msa[2].as_bytes() == b"AATGCC-CGAT");
        assert!(msa[3].as_bytes() == b"AACGCC-CGTC");
        assert!(msa[4].as_bytes() == b"AGTG-CTCGTT");
        assert!(msa[5].as_bytes() == b"AATG-CTCGTT");
        assert!(msa[6].as_bytes() == b"AATGCC-CGTT");
    }
}
