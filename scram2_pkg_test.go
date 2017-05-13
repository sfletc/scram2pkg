package scram2pkg

import (
	"github.com/montanaflynn/stats"
	"math"
	"reflect"
	"testing"
	//"fmt"

	"fmt"
)

func TestSeqLoad_single(t *testing.T) {
	var seq_files []string
	seq_files = append(seq_files, "./test_data/test_seq_1.fa")
	test_seq := SeqLoad(seq_files, "cfa", "nil", 18, 32, 2.0)
	should_be := make(map[string]*mean_se)
	var single_mean_se *mean_se
	single_mean_se = &mean_se{500000.0, 0.0}
	should_be["AAAAAAAAAAAAAAAAAAAAAAAA"] = single_mean_se
	single_mean_se = &mean_se{250000.0, 0.0}
	should_be["GGGGGGGGGGGGGGGGGGGGGGGG"] = single_mean_se
	single_mean_se = &mean_se{250000.0, 0.0}
	should_be["GGGGGGGGGGGGGGGGGGGGGGGC"] = single_mean_se
	for read, mean_ses := range test_seq {
		fmt.Println(read, mean_ses.Mean, mean_ses.Se)
	}

	eq := reflect.DeepEqual(test_seq, should_be)
	if eq == false {
		t.Error("SeqLoad not working for single seq")
	}
}

func TestSeqLoad_clean(t *testing.T) {
	var seq_files []string
	seq_files = append(seq_files, "./test_data/test_seq_5.fa")
	test_seq := SeqLoad(seq_files, "clean", "nil", 18, 32, 2.0)
	fmt.Println(test_seq )
	should_be := make(map[string]*mean_se)
	var single_mean_se *mean_se
	single_mean_se = &mean_se{500000.0, 0.0}
	should_be["AAAAAAAAAAAAAAAAAAAAAAAA"] = single_mean_se
	single_mean_se = &mean_se{250000.0, 0.0}
	should_be["GGGGGGGGGGGGGGGGGGGGGGGG"] = single_mean_se
	single_mean_se = &mean_se{250000.0, 0.0}
	should_be["GGGGGGGGGGGGGGGGGGGGGGGC"] = single_mean_se
	for read, mean_ses := range test_seq {
		fmt.Println(read, mean_ses.Mean, mean_ses.Se)
	}

	eq := reflect.DeepEqual(test_seq, should_be)
	if eq == false {
		t.Error("SeqLoad not working for single seq")
	}
}

func TestSeqLoad_fasta(t *testing.T) {
	var seq_files []string
	seq_files = append(seq_files, "./test_data/test_fasta.fasta")
	test_seq := SeqLoad(seq_files, "fa", "nil", 18, 32, 2.0)
	should_be := make(map[string]*mean_se)
	var single_mean_se *mean_se
	single_mean_se = &mean_se{500000.0, 0.0}
	should_be["AAAAAAAAAAAAAAAAAAAAAAAAA"] = single_mean_se
	single_mean_se = &mean_se{500000.0, 0.0}
	should_be["TAAAAAAAAAAAAAAAAAAAAAAAA"] = single_mean_se
	eq := reflect.DeepEqual(test_seq, should_be)
	if eq == false {
		t.Error("SeqLoad not working for single fasta file")
	}

}

//func TestSeqLoad_fa_fq(t *testing.T) {
//	var seq_files_1 []string
//	var seq_files_2 []string
//	var seq_files_3 []string
//	seq_files_1 = append(seq_files_1, "c:/Users/steve/Desktop/1RDB.fq.gz")
//	test_seq_1 := SeqLoad(seq_files_1, "fq", "nil", 18, 32, 2.0)
//	seq_files_2 = append(seq_files_2, "c:/Users/steve/Desktop/1RDB.fasta")
//	test_seq_2 := SeqLoad(seq_files_2, "fa", "nil", 18, 32, 2.0)
//	seq_files_3 = append(seq_files_3, "c:/Users/steve/Desktop/1RDB.fa")
//	test_seq_3 := SeqLoad(seq_files_3, "cfa", "nil", 18, 32, 2.0)
//	eq_1 := reflect.DeepEqual(test_seq_1, test_seq_2)
//	if eq_1 == false {
//		t.Error("fq and fasta do not match")
//	}
//	eq_2 := reflect.DeepEqual(test_seq_1, test_seq_3)
//	if eq_2 == false {
//		t.Error("cfa and fq/fasta do not match")
//	}
//}

func TestSeqLoad_multi(t *testing.T) {
	var seq_files []string
	seq_files = append(seq_files, "./test_data/test_seq_1.fa", "./test_data/test_seq_2.fa")
	test_seq := SeqLoad(seq_files, "cfa", "nil", 18, 32, 1.0)
	should_be := make(map[string]*mean_se)

	var single_mean_se *mean_se

	counts_1 := []float64{500000.0, 250000.0}
	se_1, _ := stats.StandardDeviationSample(counts_1)
	single_mean_se = &mean_se{375000.0, se_1 / math.Sqrt(2.0)}
	should_be["AAAAAAAAAAAAAAAAAAAAAAAA"] = single_mean_se

	counts_2 := []float64{250000.0, 500000.0}
	se_2, _ := stats.StandardDeviationSample(counts_2)
	single_mean_se = &mean_se{375000.0, se_2 / math.Sqrt(2.0)}
	should_be["GGGGGGGGGGGGGGGGGGGGGGGG"] = single_mean_se

	counts_3 := []float64{250000.0, 250000.0}
	se_3, _ := stats.StandardDeviationSample(counts_3)
	single_mean_se = &mean_se{250000.0, se_3 / math.Sqrt(2.0)}
	should_be["GGGGGGGGGGGGGGGGGGGGGGGC"] = single_mean_se

	eq := reflect.DeepEqual(test_seq, should_be)
	if eq == false {
		t.Error("SeqLoad not working for multi seq")
	}

}


func TestSeqLoad_multi_minCount(t *testing.T) {
	var seq_files []string
	seq_files = append(seq_files, "./test_data/test_seq_1.fa", "./test_data/test_seq_4.fa")
	test_seq := SeqLoad(seq_files, "cfa", "nil", 18, 32, 26)
	should_be := make(map[string]*mean_se)

	var single_mean_se *mean_se

	counts_1 := []float64{1000000.0, 500000.0}
	se_1, _ := stats.StandardDeviationSample(counts_1)
	single_mean_se = &mean_se{750000.0, se_1 / math.Sqrt(2.0)}
	should_be["AAAAAAAAAAAAAAAAAAAAAAAA"] = single_mean_se
	eq := reflect.DeepEqual(test_seq, should_be)
	if eq == false {
		t.Error("SeqLoad not working for multi seq")
	}

}



//func TestSeqLoad_bad_format(t *testing.T) {
//	var seq_files []string
//	seq_files = append(seq_files,"./test_data/test_seq_3_badformat.fa")
//	test_seq := SeqLoad(seq_files,18,32,1.0)
//	fmt.Println(test_seq)

//}

func TestRefLoad(t *testing.T) {
	test_ref := RefLoad("./test_data/test_ref.fa")
	var should_be []*header_ref
	ref1 := &header_ref{"ref_1", "AAAAAAAAAAAAAAAAAAAAAAAAA"}
	ref2 := &header_ref{"ref_2", "GGGGGGGGGGGGGGGGGGGGGGGGTAAAAAAAAAAAAAAAAAAAAAAAAG"}
	ref3 := &header_ref{"ref_3", ""}
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
	test_seq := SeqLoad(seq_files, "cfa", "nil", 18, 32, 1.0)
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
