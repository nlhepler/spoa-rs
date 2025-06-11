use spoa_sys::ffi;

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
    pub fn align(&mut self, sequence: &[u8], graph: &Graph) -> Alignment {
        let sequence_len = u32::try_from(sequence.len()).unwrap();
        Alignment(unsafe {
            ffi::align(
                self.0.pin_mut(),
                sequence.as_ptr() as *const i8,
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

    /// Add an `Alignment` to a `Graph`, along with its `sequence` and uniform `weight`.
    ///
    /// The `alignment` should be derived from calling `AlignmentEngine::align` with the
    /// `sequence`.
    pub fn add_alignment(&mut self, alignment: &Alignment, sequence: &[u8], weight: u32) {
        let sequence_len = u32::try_from(sequence.len()).unwrap();
        unsafe {
            ffi::add_alignment(
                self.0.pin_mut(),
                alignment.0.as_ref().unwrap(),
                sequence.as_ptr() as *const i8,
                sequence_len,
                weight,
            )
        }
    }

    /// Add an `Alignment` to a `Graph`, along with its `sequence` and per-base FASTQ `quality` scores.
    ///
    /// The `alignment` should be derived from calling `AlignmentEngine::align` with the
    /// `sequence`.
    pub fn add_alignment_with_qual(
        &mut self,
        alignment: &Alignment,
        sequence: &[u8],
        quality: &[u8],
    ) {
        let sequence_len = u32::try_from(sequence.len()).unwrap();
        let quality_len = u32::try_from(quality.len()).unwrap();
        assert!(sequence_len == quality_len);
        unsafe {
            ffi::add_alignment_with_qual(
                self.0.pin_mut(),
                alignment.0.as_ref().unwrap(),
                sequence.as_ptr() as *const i8,
                sequence_len,
                quality.as_ptr() as *const i8,
                quality_len,
            )
        }
    }

    /// Generate a consenus sequence from the partial order `Graph`.
    pub fn consensus(&mut self) -> Vec<u8> {
        Vec::from(
            ffi::generate_consensus(self.0.pin_mut())
                .as_ref()
                .unwrap()
                .as_bytes(),
        )
    }

    /// Generate a multiple sequence alignment for all sequences added to the `Graph`.
    ///
    /// If `include_consensus` is provided, the consensus sequence is provided, also aligned, at
    /// the end.
    pub fn multiple_sequence_alignment(&mut self, include_consensus: bool) -> Vec<Vec<u8>> {
        let msa = ffi::generate_multiple_sequence_alignment(self.0.pin_mut(), include_consensus);
        let mut result = vec![];
        for aln in msa.as_ref().unwrap().iter() {
            result.push(Vec::from(aln.as_bytes()));
        }
        result
    }

    /// Clear the graph
    pub fn clear(&mut self) {
        ffi::graph_clear(self.0.pin_mut())
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
    const SMALL_SEQS: &[&[u8]] = &[
        b"ATTGCCCGTT",
        b"AATGCCGTT",
        b"AATGCCCGAT",
        b"AACGCCCGTC",
        b"AGTGCTCGTT",
        b"AATGCTCGTT",
    ];

    #[test]
    fn consensus_small() {
        let mut eng = AlignmentEngine::new(AlignmentType::kNW, 5, -4, -3, -1, -3, -1);
        let mut graph = Graph::new();

        for seq in SMALL_SEQS {
            let aln = eng.align(seq, &graph);
            graph.add_alignment(&aln, seq, 1);
        }

        let consensus = graph.consensus();
        assert!(consensus.as_slice() == b"AATGCCCGTT");
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
            let aln = eng.align(record.seq.as_bytes(), &graph);
            graph.add_alignment_with_qual(&aln, record.seq.as_bytes(), record.qual.as_bytes());
        }

        let consensus = graph.consensus();
        let expected = b"AATGATGCGCTTTGTTGGCGCGGTGGCTTGATGCAGGGGCTAATCGACCTCTGGCAACCACTTTTCCATGACAGGAGTTGAATATGGCATTCAGTAATCCCTTCGATGATCCGCAGGGAGCGTTTTACATATTGCGCAATGCGCAGGGGCAATTCAGTCTGTGGCCGCAACAATGCGTCTTACCGGCAGGCTGGGACATTGTGTGTCAGCCGCAGTCACAGGCGTCCTGCCAGCAGTGGCTGGAAGCCCACTGGCGTACTCTGACACCGACGAATTTTACCCAGTTGCAGGAGGCACAATGAGCCAGCATTTACCTTTGGTCGCCGCACAGCCCGGCATCTGGATGGCAGAAAAACTGTCAGAATTACCCTCCGCCTGGAGCGTGGCGCATTACGTTGAGTTAACCGGAGAGGTTGATTCGCCATTACTGGCCCGCGCGGTGGTTGCCGGACTAGCGCAAGCAGATACGCTTTACACGCGCAACCAAGGATTTCGG";
        assert!(consensus.as_slice() == expected);
    }

    #[test]
    fn msa_small() {
        let mut eng = AlignmentEngine::new(AlignmentType::kNW, 5, -4, -3, -1, -3, -1);
        let mut graph = Graph::new();

        for seq in SMALL_SEQS {
            let aln = eng.align(seq, &graph);
            graph.add_alignment(&aln, seq, 1);
        }

        let msa = graph.multiple_sequence_alignment(true);
        assert!(msa.len() == 7);
        assert!(msa[0].as_slice() == b"ATTGCC-CGTT");
        assert!(msa[1].as_slice() == b"AATG-C-CGTT");
        assert!(msa[2].as_slice() == b"AATGCC-CGAT");
        assert!(msa[3].as_slice() == b"AACGCC-CGTC");
        assert!(msa[4].as_slice() == b"AGTG-CTCGTT");
        assert!(msa[5].as_slice() == b"AATG-CTCGTT");
        assert!(msa[6].as_slice() == b"AATGCC-CGTT");

        graph.clear();

        for seq in SMALL_SEQS.iter().rev() {
            let aln = eng.align(seq, &graph);
            graph.add_alignment(&aln, seq, 1);
        }

        let msa = graph.multiple_sequence_alignment(true);
        assert!(msa.len() == 7);
        assert!(msa[0].as_slice() == b"AATGCTCGTT");
        assert!(msa[1].as_slice() == b"AGTGCTCGTT");
        assert!(msa[2].as_slice() == b"AACGCCCGTC");
        assert!(msa[3].as_slice() == b"AATGCCCGAT");
        assert!(msa[4].as_slice() == b"AATGC-CGTT");
        assert!(msa[5].as_slice() == b"ATTGCCCGTT");
        assert!(msa[6].as_slice() == b"AATGCCCGTT");
    }
}
