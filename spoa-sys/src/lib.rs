#[cxx::bridge(namespace = "spoa")]
pub mod ffi {
    /// Enumerates the possible alignment types
    #[repr(u32)]
    enum AlignmentType {
        /// Smith-Waterman alignment
        kSW,
        /// Needleman-Wunsch alignment
        kNW,
        /// Overlap alignment
        kOV,
    }

    unsafe extern "C++" {
        include!("spoa-sys/include/bindings.hpp");

        type Alignment;
        type AlignmentEngine;
        type AlignmentType;
        type Graph;

        fn create_graph() -> UniquePtr<Graph>;
        /// # Safety
        /// this function is unsafe because cxx says pointer arguments are unsafe
        unsafe fn add_alignment_with_qual(
            graph: Pin<&mut Graph>,
            alignment: &Alignment,
            sequence: *const c_char,
            sequence_len: u32,
            quality: *const c_char,
            quality_len: u32,
        );
        unsafe fn add_alignment(
            graph: Pin<&mut Graph>,
            alignment: &Alignment,
            sequence: *const c_char,
            sequence_len: u32,
            weight: u32,
        );
        fn generate_consensus(graph: Pin<&mut Graph>) -> UniquePtr<CxxString>;
        fn generate_consensus_with_min_coverage(
            graph: Pin<&mut Graph>,
            min_coverage: i32,
        ) -> UniquePtr<CxxString>;
        fn generate_multiple_sequence_alignment(
            graph: Pin<&mut Graph>,
            include_consensus: bool,
        ) -> UniquePtr<CxxVector<CxxString>>;
        fn graph_clear(graph: Pin<&mut Graph>);

        fn create_alignment_engine(
            typ: AlignmentType,
            m: i8,
            n: i8,
            g: i8,
            e: i8,
            q: i8,
            c: i8,
        ) -> UniquePtr<AlignmentEngine>;
        /// # Safety
        /// this function is unsafe because cxx says pointer arguments are unsafe
        unsafe fn align(
            alignment_engine: Pin<&mut AlignmentEngine>,
            sequence: *const c_char,
            sequence_len: u32,
            graph: &Graph,
        ) -> UniquePtr<Alignment>;
    }
}
