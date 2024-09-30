use timsseek::fragment_mass::elution_group_converter::SequenceToElutionGroupConverter;
use timsseek::fragment_mass::fragment_mass_builder::FragmentMassBuilder;

// timseek --fasta asdasdad --config asdad.json --out_dir outputs # should generate a
// 'results.sqlite' with 1 score per unique peptide+charge in the fasta file.
// "Score" is used loosely here, its the combination of score + RT

fn main() {
    let seq = "PEPTIDEPINK/2";
    // let fragment_mass_builder = FragmentMassBuilder::default();
    // let fragment_mzs = fragment_mass_builder
    //     .fragment_mzs_from_proforma(seq)
    //     .unwrap();

    // println!("{:?}", fragment_mzs);
    // println!("Hello, world!");

    let def_converter = SequenceToElutionGroupConverter::default();
    let out = def_converter.convert_sequences(&[seq]);
    println!("{:?}", out);
}
