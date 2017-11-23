package scramPkg

import (
	"fmt"
	"github.com/montanaflynn/stats"
	"math"
	"reflect"
	"testing"
)

func TestSeqLoad_single(t *testing.T) {
	var seq_files []string
	seq_files = append(seq_files, "./test_data/test_seq_1.fa")
	test_seq := SeqLoad(seq_files, "cfa", "nil", 18, 32, 2.0, false)
	should_be := make(map[string]interface{})
	var single_mean_se *meanSe
	single_mean_se = &meanSe{500000.0, 0.0}
	should_be["AAAAAAAAAAAAAAAAAAAAAAAA"] = single_mean_se
	single_mean_se = &meanSe{250000.0, 0.0}
	should_be["GGGGGGGGGGGGGGGGGGGGGGGG"] = single_mean_se
	single_mean_se = &meanSe{250000.0, 0.0}
	should_be["GGGGGGGGGGGGGGGGGGGGGGGC"] = single_mean_se
	fmt.Println(test_seq)
	for read, mean_ses := range test_seq {
		fmt.Println(read, mean_ses.(*meanSe).Mean, mean_ses.(*meanSe).Se)
	}
	eq := reflect.DeepEqual(test_seq, should_be)
	if eq == false {
		t.Error("SeqLoad not working for single seq")
	}
}

func TestIndvSeqLoad_single(t *testing.T) {

	var seq_files []string
	seq_files = append(seq_files, "./test_data/test_seq_1.fa")
	test_seq, _ := IndvSeqLoad(seq_files, "cfa", "nil", 18, 32, 2.0, false)
	should_be := make(map[string]interface{})
	var indv_counts *[]float64
	indv_counts = &[]float64{500000.0}
	should_be["AAAAAAAAAAAAAAAAAAAAAAAA"] = indv_counts
	indv_counts = &[]float64{250000.0}
	should_be["GGGGGGGGGGGGGGGGGGGGGGGG"] = indv_counts
	indv_counts = &[]float64{250000.0}
	should_be["GGGGGGGGGGGGGGGGGGGGGGGC"] = indv_counts
	eq := reflect.DeepEqual(test_seq, should_be)
	if eq == false {
		t.Error("SeqLoad not working for single seq")
	}
}

func TestSeqLoad_clean(t *testing.T) {
	var seq_files []string
	seq_files = append(seq_files, "./test_data/test_seq_5.fa")
	test_seq := SeqLoad(seq_files, "clean", "nil", 18, 32, 2.0, false)
	fmt.Println(test_seq)
	should_be := make(map[string]interface{})
	var single_mean_se *meanSe
	single_mean_se = &meanSe{500000.0, 0.0}
	should_be["AAAAAAAAAAAAAAAAAAAAAAAA"] = single_mean_se
	single_mean_se = &meanSe{250000.0, 0.0}
	should_be["GGGGGGGGGGGGGGGGGGGGGGGG"] = single_mean_se
	single_mean_se = &meanSe{250000.0, 0.0}
	should_be["GGGGGGGGGGGGGGGGGGGGGGGC"] = single_mean_se
	for read, mean_ses := range test_seq {
		fmt.Println(read, mean_ses.(*meanSe).Mean, mean_ses.(*meanSe).Se)
	}
	eq := reflect.DeepEqual(test_seq, should_be)
	if eq == false {
		t.Error("SeqLoad not working for single seq")
	}
}

func TestSeqLoad_fasta(t *testing.T) {
	var seq_files []string
	seq_files = append(seq_files, "./test_data/test_fasta.fasta")
	test_seq := SeqLoad(seq_files, "fa", "nil", 18, 32, 2.0, false)
	should_be := make(map[string]interface{})
	var single_mean_se *meanSe
	single_mean_se = &meanSe{500000.0, 0.0}
	should_be["AAAAAAAAAAAAAAAAAAAAAAAAA"] = single_mean_se
	single_mean_se = &meanSe{500000.0, 0.0}
	should_be["TAAAAAAAAAAAAAAAAAAAAAAAA"] = single_mean_se
	eq := reflect.DeepEqual(test_seq, should_be)
	if eq == false {
		t.Error("SeqLoad not working for single fasta file")
	}

}

func TestSeqLoad_multi(t *testing.T) {
	var seq_files []string
	seq_files = append(seq_files, "./test_data/test_seq_1.fa", "./test_data/test_seq_2.fa")
	test_seq := SeqLoad(seq_files, "cfa", "nil", 18, 32, 1.0, false)
	fmt.Println(test_seq)
	fmt.Println(test_seq["ACGCTGATGCATGCATCGACTAGC"])
	should_be := make(map[string]interface{})

	var single_mean_se *meanSe

	counts_1 := []float64{500000.0, 250000.0}
	se_1, _ := stats.StandardDeviationSample(counts_1)
	single_mean_se = &meanSe{375000.0, se_1 / math.Sqrt(2.0)}
	should_be["AAAAAAAAAAAAAAAAAAAAAAAA"] = single_mean_se

	counts_2 := []float64{250000.0, 500000.0}
	se_2, _ := stats.StandardDeviationSample(counts_2)
	single_mean_se = &meanSe{375000.0, se_2 / math.Sqrt(2.0)}
	should_be["GGGGGGGGGGGGGGGGGGGGGGGG"] = single_mean_se

	counts_3 := []float64{250000.0, 250000.0}
	se_3, _ := stats.StandardDeviationSample(counts_3)
	single_mean_se = &meanSe{250000.0, se_3 / math.Sqrt(2.0)}
	should_be["GGGGGGGGGGGGGGGGGGGGGGGC"] = single_mean_se

	eq := reflect.DeepEqual(test_seq, should_be)
	if eq == false {
		t.Error("SeqLoad not working for multi seq")
	}

}

func TestIndvSeqLoad_multi(t *testing.T) {
	var seq_files []string
	seq_files = append(seq_files, "./test_data/test_seq_1.fa", "./test_data/test_seq_2.fa")
	test_seq, _ := IndvSeqLoad(seq_files, "cfa", "nil", 18, 32, 1.0, false)
	fmt.Println(test_seq)
	fmt.Println(test_seq["ACGCTGATGCATGCATCGACTAGC"])
	should_be := make(map[string]interface{})

	var indv_counts *[]float64

	indv_counts = &[]float64{500000.0, 250000.0}
	should_be["AAAAAAAAAAAAAAAAAAAAAAAA"] = indv_counts

	indv_counts = &[]float64{250000.0, 500000.0}
	should_be["GGGGGGGGGGGGGGGGGGGGGGGG"] = indv_counts

	indv_counts = &[]float64{250000.0, 250000.0}
	should_be["GGGGGGGGGGGGGGGGGGGGGGGC"] = indv_counts

	for read, counts := range test_seq {
		fmt.Println(read, counts)
	}

	eq := reflect.DeepEqual(test_seq, should_be)
	if eq == false {
		t.Error("SeqLoad not working for multi seq")
	}

}

func TestSeqLoad_multi_minCount(t *testing.T) {
	var seq_files []string
	seq_files = append(seq_files, "./test_data/test_seq_1.fa", "./test_data/test_seq_4.fa")
	test_seq := SeqLoad(seq_files, "cfa", "nil", 18, 32, 26, false)
	fmt.Println(test_seq)
	should_be := make(map[string]interface{})

	var single_mean_se *meanSe

	counts_1 := []float64{1000000.0, 500000.0}
	se_1, _ := stats.StandardDeviationSample(counts_1)
	single_mean_se = &meanSe{750000.0, se_1 / math.Sqrt(2.0)}
	should_be["AAAAAAAAAAAAAAAAAAAAAAAA"] = single_mean_se
	eq := reflect.DeepEqual(test_seq, should_be)
	if eq == false {
		t.Error("SeqLoad not working for multi seq")
	}

}

func TestRefLoad(t *testing.T) {
	test_ref := RefLoad("./test_data/test_ref.fa")
	var should_be []*headerRef
	ref1 := &headerRef{"ref_1", "AAAAAAAAAAAAAAAAAAAAAAAAA", "TTTTTTTTTTTTTTTTTTTTTTTTT"}
	ref2 := &headerRef{"ref_2", "GGGGGGGGGGGGGGGGGGGGGGGGTAAAAAAAAAAAAAAAAAAAAAAAAG", "CTTTTTTTTTTTTTTTTTTTTTTTTACCCCCCCCCCCCCCCCCCCCCCCC"}
	ref3 := &headerRef{"ref_3", "", ""}
	should_be = append(should_be, ref1, ref2, ref3)
	fmt.Println(test_ref)
	if len(test_ref) != len(should_be) {
		t.Error("Wrong no of refs in test_ref.fa")
	}
	for i := 0; i < len(test_ref); i++ {
		a := test_ref[i].header
		b := should_be[i].header
		if a != b {
			t.Error("Headers dont't match")
		}
		if test_ref[i].seq != should_be[i].seq {
			fmt.Println(test_ref[i].seq)
			fmt.Println(should_be[i].seq)
			t.Error("Seqs dont't match")
		}
	}
}

func TestAlign(t *testing.T) {
	test_ref := RefLoad("./test_data/test_ref_align.fa")

	var seq_files []string
	seq_files = append(seq_files, "./test_data/test_seq_1.fa")
	test_seq := SeqLoad(seq_files, "cfa", "nil", 18, 32, 1.0, false)
	test_align := AlignReads(test_seq, test_ref, 24)
	pos_1 := []int{1, 2}
	single_align_1 := map[string][]int{"AAAAAAAAAAAAAAAAAAAAAAAA": pos_1}
	pos_2 := []int{1}
	pos_3 := []int{26}
	single_align_2 := map[string][]int{"GGGGGGGGGGGGGGGGGGGGGGGG": pos_2, "AAAAAAAAAAAAAAAAAAAAAAAA": pos_3}
	pos_4 := []int{-2, -1}
	single_align_3 := map[string][]int{"AAAAAAAAAAAAAAAAAAAAAAAA": pos_4}
	should_be := map[string]map[string][]int{"ref_1": single_align_1, "ref_2": single_align_2, "ref_3": single_align_3}
	eq := reflect.DeepEqual(test_align, should_be)
	if eq == false {
		t.Error("Alignments are not equal")
	}

}

func TestCompareAlign(t *testing.T) {
	test_ref := RefLoad("./test_data/test_ref_align.fa")

	var seq_files []string
	seq_files = append(seq_files, "./test_data/test_seq_1.fa")
	test_seq := SeqLoad(seq_files, "cfa", "nil", 18, 32, 1.0, false)
	test_align_1 := AlignReads(test_seq, test_ref, 24)
	test_align_2 := AlignReads(test_seq, test_ref, 24)
	fmt.Println(test_align_1, test_align_2)
	test_align_1_split := CompareSplitCounts(test_align_1, test_seq)
	test_align_2_split := CompareSplitCounts(test_align_2, test_seq)
	fmt.Println(test_align_1_split, test_align_2_split)
	test_align_1_nosplit := CompareNoSplitCounts(test_align_1, test_seq)
	test_align_2_nosplit := CompareNoSplitCounts(test_align_2, test_seq)
	fmt.Println(test_align_1_nosplit, test_align_2_nosplit)
	compare_split := Compare(test_align_1_split, test_align_2_split)
	compare_no_split := Compare(test_align_1_nosplit, test_align_2_nosplit)
	fmt.Println(compare_split)

	compare_split_should_be := make(map[string]interface{})
	compare_split_should_be["ref_1"] = compMeanSeOutput{}
	compare_split_should_be["ref_1"] = compare_split_should_be["ref_1"].(compMeanSeOutput).append(200000.0, 0.0, 200000.0, 0.0)
	compare_split_should_be["ref_2"] = compMeanSeOutput{}
	compare_split_should_be["ref_2"] = compare_split_should_be["ref_2"].(compMeanSeOutput).append(350000.0, 0.0, 350000.0, 0.0)
	compare_split_should_be["ref_3"] = compMeanSeOutput{}
	compare_split_should_be["ref_3"] = compare_split_should_be["ref_3"].(compMeanSeOutput).append(200000.0, 0.0, 200000.0, 0.0)

	fmt.Println(compare_split_should_be)
	eq := reflect.DeepEqual(compare_split, compare_split_should_be)
	if eq == false {
		t.Error("Compare split alignments not equal")
	}

	compare_no_split_should_be := make(map[string]interface{})
	compare_no_split_should_be["ref_1"] = compMeanSeOutput{}
	compare_no_split_should_be["ref_1"] = compare_no_split_should_be["ref_1"].(compMeanSeOutput).append(1000000.0, 0.0, 1000000.0, 0.0)
	compare_no_split_should_be["ref_2"] = compMeanSeOutput{}
	compare_no_split_should_be["ref_2"] = compare_no_split_should_be["ref_2"].(compMeanSeOutput).append(750000.0, 0.0, 750000.0, 0.0)
	compare_no_split_should_be["ref_3"] = compMeanSeOutput{}
	compare_no_split_should_be["ref_3"] = compare_no_split_should_be["ref_3"].(compMeanSeOutput).append(1000000.0, 0.0, 1000000.0, 0.0)
	fmt.Println(compare_no_split)
	fmt.Println(compare_no_split_should_be)
	eq2 := reflect.DeepEqual(compare_no_split, compare_no_split_should_be)
	if eq2 == false {
		t.Error("Compare no split alignments not equal")
	}
}

func TestIndvCompareAlign(t *testing.T) {
	test_ref := RefLoad("./test_data/test_ref_align.fa")

	var seq_files []string
	seq_files = append(seq_files, "./test_data/test_seq_1.fa", "./test_data/test_seq_1.fa")
	test_seq, _ := IndvSeqLoad(seq_files, "cfa", "nil", 18, 32, 1.0, false)
	for key, value := range test_seq {
		fmt.Println(key, value)
	}
	test_align_1 := AlignReads(test_seq, test_ref, 24)
	test_align_2 := AlignReads(test_seq, test_ref, 24)
	fmt.Println(test_align_1, test_align_2)
	test_align_1_split := CompareSplitCounts(test_align_1, test_seq)
	test_align_2_split := CompareSplitCounts(test_align_2, test_seq)
	fmt.Println(test_align_1_split, test_align_2_split)
	test_align_1_nosplit := CompareNoSplitCounts(test_align_1, test_seq)
	test_align_2_nosplit := CompareNoSplitCounts(test_align_2, test_seq)
	fmt.Println(test_align_1_nosplit, test_align_2_nosplit)
	compare_split := Compare(test_align_1_split, test_align_2_split)
	compare_no_split := Compare(test_align_1_nosplit, test_align_2_nosplit)
	fmt.Println(compare_split)

	compare_split_should_be := make(map[string]interface{})
	compare_split_should_be["ref_1"] = countsOutput{}
	compare_split_should_be["ref_1"] = compare_split_should_be["ref_1"].(countsOutput).append(200000.0, 200000.0, 200000.0, 200000.0)
	compare_split_should_be["ref_2"] = countsOutput{}
	compare_split_should_be["ref_2"] = compare_split_should_be["ref_2"].(countsOutput).append(350000.0, 350000.0, 350000.0, 350000.0)
	compare_split_should_be["ref_3"] = countsOutput{}
	compare_split_should_be["ref_3"] = compare_split_should_be["ref_3"].(countsOutput).append(200000.0, 200000.0, 200000.0, 200000.0)

	fmt.Println(compare_split_should_be)
	eq := reflect.DeepEqual(compare_split, compare_split_should_be)
	if eq == false {
		t.Error("Indv Compare split alignments not equal")
	}
	fmt.Println(compare_no_split)
	compare_no_split_should_be := make(map[string]interface{})
	compare_no_split_should_be["ref_1"] = countsOutput{}
	compare_no_split_should_be["ref_1"] = compare_no_split_should_be["ref_1"].(countsOutput).append(1000000.0, 1000000.0, 1000000.0, 1000000.0)
	compare_no_split_should_be["ref_2"] = countsOutput{}
	compare_no_split_should_be["ref_2"] = compare_no_split_should_be["ref_2"].(countsOutput).append(750000.0, 750000.0, 750000.0, 750000.0)
	compare_no_split_should_be["ref_3"] = countsOutput{}
	compare_no_split_should_be["ref_3"] = compare_no_split_should_be["ref_3"].(countsOutput).append(1000000.0, 1000000.0, 1000000.0, 1000000.0)
	fmt.Println(compare_no_split)
	fmt.Println(compare_no_split_should_be)
	eq2 := reflect.DeepEqual(compare_no_split, compare_no_split_should_be)
	if eq2 == false {
		t.Error("Indv Compare split alignments not equal")
	}
}

func TestProfileAlign(t *testing.T) {
	test_ref := RefLoad("./test_data/test_ref_align.fa")

	var seq_files []string
	seq_files = append(seq_files, "./test_data/test_seq_1.fa")
	test_seq := SeqLoad(seq_files, "cfa", "nil", 18, 32, 1.0, false)
	test_align_1 := AlignReads(test_seq, test_ref, 24)

	fmt.Println(test_align_1)
	test_align_1_split := ProfileSplit(test_align_1, test_seq)
	for ref, value := range test_align_1_split {
		fmt.Println(ref)
		for _, i := range *value.(*singleAlignments) {
			fmt.Println(i, i.Alignments)
		}
	}
	should_be_split := make(map[string]interface{})
	combinedAlignments_1 := singleAlignments{}
	alignment_1 := singleAlignment{"AAAAAAAAAAAAAAAAAAAAAAAA", 5, 1, "+", &meanSe{100000, 0}}
	alignment_2 := singleAlignment{"AAAAAAAAAAAAAAAAAAAAAAAA", 5, 2, "+", &meanSe{100000, 0}}
	combinedAlignments_1 = append(combinedAlignments_1, &alignment_1, &alignment_2)
	should_be_split["ref_1"] = &combinedAlignments_1

	combinedAlignments_2 := singleAlignments{}
	alignment_3 := singleAlignment{"GGGGGGGGGGGGGGGGGGGGGGGG", 1, 1, "+", &meanSe{250000, 0}}
	alignment_4 := singleAlignment{"AAAAAAAAAAAAAAAAAAAAAAAA", 5, 26, "+", &meanSe{100000, 0}}
	combinedAlignments_2 = append(combinedAlignments_2, &alignment_3, &alignment_4)
	should_be_split["ref_2"] = &combinedAlignments_2

	combinedAlignments_3 := singleAlignments{}
	alignment_5 := singleAlignment{"AAAAAAAAAAAAAAAAAAAAAAAA", 5, 1, "-", &meanSe{100000, 0}}
	alignment_6 := singleAlignment{"AAAAAAAAAAAAAAAAAAAAAAAA", 5, 2, "-", &meanSe{100000, 0}}
	combinedAlignments_3 = append(combinedAlignments_3, &alignment_5, &alignment_6)
	should_be_split["ref_3"] = &combinedAlignments_3
	for ref, value := range should_be_split {
		fmt.Println(ref)
		for _, i := range *value.(*singleAlignments) {
			fmt.Println(i, i.Alignments)
		}
	}

	eq1 := reflect.DeepEqual(should_be_split, test_align_1_split)
	if eq1 == false {
		t.Error("Profile alignment is incorrect")
	}

	test_align_1_no_split := ProfileNoSplit(test_align_1, test_seq)

	should_be_no_split := make(map[string]interface{})
	combinedAlignments_4 := singleAlignments{}
	alignment_7 := singleAlignment{"AAAAAAAAAAAAAAAAAAAAAAAA", 5, 1, "+", &meanSe{500000, 0}}
	alignment_8 := singleAlignment{"AAAAAAAAAAAAAAAAAAAAAAAA", 5, 2, "+", &meanSe{500000, 0}}
	combinedAlignments_4 = append(combinedAlignments_4, &alignment_7, &alignment_8)
	should_be_no_split["ref_1"] = &combinedAlignments_4

	combinedAlignments_5 := singleAlignments{}
	alignment_9 := singleAlignment{"GGGGGGGGGGGGGGGGGGGGGGGG", 1, 1, "+", &meanSe{250000, 0}}
	alignment_10 := singleAlignment{"AAAAAAAAAAAAAAAAAAAAAAAA", 5, 26, "+", &meanSe{500000, 0}}
	combinedAlignments_5 = append(combinedAlignments_5, &alignment_9, &alignment_10)
	should_be_no_split["ref_2"] = &combinedAlignments_5

	combinedAlignments_6 := singleAlignments{}
	alignment_11 := singleAlignment{"AAAAAAAAAAAAAAAAAAAAAAAA", 5, 1, "-", &meanSe{500000, 0}}
	alignment_12 := singleAlignment{"AAAAAAAAAAAAAAAAAAAAAAAA", 5, 2, "-", &meanSe{500000, 0}}
	combinedAlignments_6 = append(combinedAlignments_6, &alignment_11, &alignment_12)
	should_be_no_split["ref_3"] = &combinedAlignments_6

	eq2 := reflect.DeepEqual(should_be_no_split, test_align_1_no_split)
	if eq2 == false {
		t.Error("Profile no split alignment is incorrect")
	}
}

func TestProfileAlignIndv(t *testing.T) {
	test_ref := RefLoad("./test_data/test_ref_align.fa")

	var seq_files []string
	seq_files = append(seq_files, "./test_data/test_seq_1.fa")
	test_seq, _ := IndvSeqLoad(seq_files, "cfa", "nil", 18, 32, 1.0, false)
	test_align_1 := AlignReads(test_seq, test_ref, 24)

	fmt.Println(test_align_1)
	test_align_1_split := ProfileSplit(test_align_1, test_seq)
	for ref, value := range test_align_1_split {
		fmt.Println(ref)
		for _, i := range *value.(*singleAlignments) {
			fmt.Println(i, i.Alignments)
		}
	}
	should_be_split := make(map[string]interface{})
	combinedAlignments_1 := singleAlignments{}
	alignment_1 := singleAlignment{"AAAAAAAAAAAAAAAAAAAAAAAA", 5, 1, "+", &[]float64{100000}}
	alignment_2 := singleAlignment{"AAAAAAAAAAAAAAAAAAAAAAAA", 5, 2, "+", &[]float64{100000}}
	combinedAlignments_1 = append(combinedAlignments_1, &alignment_1, &alignment_2)
	should_be_split["ref_1"] = &combinedAlignments_1

	combinedAlignments_2 := singleAlignments{}
	alignment_3 := singleAlignment{"GGGGGGGGGGGGGGGGGGGGGGGG", 1, 1, "+", &[]float64{250000}}
	alignment_4 := singleAlignment{"AAAAAAAAAAAAAAAAAAAAAAAA", 5, 26, "+", &[]float64{100000}}
	combinedAlignments_2 = append(combinedAlignments_2, &alignment_3, &alignment_4)
	should_be_split["ref_2"] = &combinedAlignments_2

	combinedAlignments_3 := singleAlignments{}
	alignment_5 := singleAlignment{"AAAAAAAAAAAAAAAAAAAAAAAA", 5, 1, "-", &[]float64{100000}}
	alignment_6 := singleAlignment{"AAAAAAAAAAAAAAAAAAAAAAAA", 5, 2, "-", &[]float64{100000}}
	combinedAlignments_3 = append(combinedAlignments_3, &alignment_5, &alignment_6)
	should_be_split["ref_3"] = &combinedAlignments_3
	for ref, value := range should_be_split {
		fmt.Println(ref)
		for _, i := range *value.(*singleAlignments) {
			fmt.Println(i, i.Alignments)
		}
	}

	eq1 := reflect.DeepEqual(should_be_split, test_align_1_split)
	if eq1 == false {
		t.Error("Profile alignment is incorrect")
	}

	test_align_1_no_split := ProfileNoSplit(test_align_1, test_seq)

	should_be_no_split := make(map[string]interface{})
	combinedAlignments_4 := singleAlignments{}
	alignment_7 := singleAlignment{"AAAAAAAAAAAAAAAAAAAAAAAA", 5, 1, "+", &[]float64{500000}}
	alignment_8 := singleAlignment{"AAAAAAAAAAAAAAAAAAAAAAAA", 5, 2, "+", &[]float64{500000}}
	combinedAlignments_4 = append(combinedAlignments_4, &alignment_7, &alignment_8)
	should_be_no_split["ref_1"] = &combinedAlignments_4

	combinedAlignments_5 := singleAlignments{}
	alignment_9 := singleAlignment{"GGGGGGGGGGGGGGGGGGGGGGGG", 1, 1, "+", &[]float64{250000}}
	alignment_10 := singleAlignment{"AAAAAAAAAAAAAAAAAAAAAAAA", 5, 26, "+", &[]float64{500000}}
	combinedAlignments_5 = append(combinedAlignments_5, &alignment_9, &alignment_10)
	should_be_no_split["ref_2"] = &combinedAlignments_5

	combinedAlignments_6 := singleAlignments{}
	alignment_11 := singleAlignment{"AAAAAAAAAAAAAAAAAAAAAAAA", 5, 1, "-", &[]float64{500000}}
	alignment_12 := singleAlignment{"AAAAAAAAAAAAAAAAAAAAAAAA", 5, 2, "-", &[]float64{500000}}
	combinedAlignments_6 = append(combinedAlignments_6, &alignment_11, &alignment_12)
	should_be_no_split["ref_3"] = &combinedAlignments_6

	eq2 := reflect.DeepEqual(should_be_no_split, test_align_1_no_split)
	if eq2 == false {
		t.Error("Profile no split alignment is incorrect")
	}
}

func TestMirAlign(t *testing.T) {
	test_mir_ref := MirLoad("./test_data/test_mir_align.fa")
	var seq_files []string
	seq_files = append(seq_files, "./test_data/test_seq_6.fa")
	test_seq := SeqLoad(seq_files, "cfa", "nil", 1, 32, 1.0, false)
	test_align_1 := AlignMirnas(test_seq, test_mir_ref)
	test_align_2 := AlignMirnas(test_seq, test_mir_ref)
	noSplitMap := MirnaCompare(test_align_1, test_align_2, true)
	splitMap := MirnaCompare(test_align_1, test_align_2, false)
	fmt.Println(noSplitMap, splitMap)
	shouldBeNoSplit := map[string]interface{}{"mir_1": compMeanSeOutput{[]float64{500000, 0, 500000, 0}},
		"mir_3": compMeanSeOutput{[]float64{500000, 0, 500000, 0}},
		"mir_2": compMeanSeOutput{[]float64{250000, 0, 250000, 0}}}
	shouldBeSplit := map[string]interface{}{"mir_1": compMeanSeOutput{[]float64{250000, 0, 250000, 0}},
		"mir_3": compMeanSeOutput{[]float64{250000, 0, 250000, 0}},
		"mir_2": compMeanSeOutput{[]float64{250000, 0, 250000, 0}}}
	fmt.Println(shouldBeNoSplit, shouldBeSplit)
	eq_1 := reflect.DeepEqual(noSplitMap, shouldBeNoSplit)
	if eq_1 == false {
		t.Error("noSplit mir alignments are not equal")
	}
	eq_2 := reflect.DeepEqual(splitMap, shouldBeSplit)
	if eq_2 == false {
		t.Error("split mir alignments are not equal")
	}
}

func TestIndvMirAlign(t *testing.T) {
	test_mir_ref := MirLoad("./test_data/test_mir_align.fa")
	var seq_files []string
	seq_files = append(seq_files, "./test_data/test_seq_6.fa")
	test_seq, _ := IndvSeqLoad(seq_files, "cfa", "nil", 1, 32, 1.0, false)
	test_align_1 := AlignMirnas(test_seq, test_mir_ref)
	test_align_2 := AlignMirnas(test_seq, test_mir_ref)
	noSplitMap := MirnaCompare(test_align_1, test_align_2, true)
	splitMap := MirnaCompare(test_align_1, test_align_2, false)
	fmt.Println(noSplitMap, splitMap)
	shouldBeNoSplit := map[string]interface{}{"mir_1": countsOutput{[]float64{500000, 500000}},
		"mir_3": countsOutput{[]float64{500000, 500000}},
		"mir_2": countsOutput{[]float64{250000, 250000}}}
	shouldBeSplit := map[string]interface{}{"mir_1": countsOutput{[]float64{250000, 250000}},
		"mir_3": countsOutput{[]float64{250000, 250000}},
		"mir_2": countsOutput{[]float64{250000, 250000}}}
	fmt.Println(shouldBeNoSplit, shouldBeSplit)
	eq_1 := reflect.DeepEqual(noSplitMap, shouldBeNoSplit)
	if eq_1 == false {
		t.Error("noSplit mir alignments are not equal")
	}
	eq_2 := reflect.DeepEqual(splitMap, shouldBeSplit)
	if eq_2 == false {
		t.Error("split mir alignments are not equal")
	}
}
