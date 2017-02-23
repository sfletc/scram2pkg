package scram2pkg

import (
	"testing"
	"reflect"
	"github.com/montanaflynn/stats"
	"math"
	//"fmt"

	"fmt"

)



func TestSeqLoad_single(t *testing.T) {
	var seq_files []string
	seq_files = append(seq_files,"./test_data/test_seq_1.fa")
	test_seq := SeqLoad(seq_files,18,32,1.0)
	should_be := make(map[string]*mean_se)
	var single_mean_se *mean_se
	single_mean_se = &mean_se{500000.0,0.0}
	should_be["AAAAAAAAAAAAAAAAAAAAAAAA"]=single_mean_se
	single_mean_se = &mean_se{250000.0,0.0}
	should_be["GGGGGGGGGGGGGGGGGGGGGGGG"]=single_mean_se
	single_mean_se = &mean_se{250000.0,0.0}
	should_be["GGGGGGGGGGGGGGGGGGGGGGGC"]=single_mean_se
	eq:=reflect.DeepEqual(test_seq,should_be)
	if eq == false {
		t.Error("SeqLoad not working for single seq")
	}

}

func TestSeqLoad_multi(t *testing.T) {
	var seq_files []string
	seq_files = append(seq_files,"./test_data/test_seq_1.fa", "./test_data/test_seq_2.fa")
	test_seq := SeqLoad(seq_files,18,32,1.0)
	should_be := make(map[string]*mean_se)

	var single_mean_se *mean_se

	counts_1 := []float64{500000.0,250000.0}
	se_1,_:=stats.StandardDeviationSample(counts_1)
	single_mean_se = &mean_se{375000.0,se_1/math.Sqrt(2.0)}
	should_be["AAAAAAAAAAAAAAAAAAAAAAAA"]=single_mean_se

	counts_2 := []float64{250000.0,500000.0}
	se_2,_:=stats.StandardDeviationSample(counts_2)
	single_mean_se = &mean_se{375000.0,se_2/math.Sqrt(2.0)}
	should_be["GGGGGGGGGGGGGGGGGGGGGGGG"]=single_mean_se

	counts_3 := []float64{250000.0,250000.0}
	se_3,_:=stats.StandardDeviationSample(counts_3)
	single_mean_se = &mean_se{250000.0,se_3/math.Sqrt(2.0)}
	should_be["GGGGGGGGGGGGGGGGGGGGGGGC"]=single_mean_se

	eq:=reflect.DeepEqual(test_seq,should_be)
	if eq == false {
		t.Error("SeqLoad not working for multi seq")
	}

}

//func TestSeqLoad_baf_format(t *testing.T) {
//	var seq_files []string
//	seq_files = append(seq_files,"./test_data/test_seq_3_badformat.fa")
//	test_seq := SeqLoad(seq_files,18,32,1.0)
//	fmt.Println(test_seq)


//}

func TestRefLoad(t *testing.T) {
	test_ref := RefLoad("./test_data/test_ref.fa")
	var should_be []*header_ref
	ref1 := &header_ref{"ref_1","AAAAAAAAAAAAAAAAAAAAAAAAA"}
	ref2 := &header_ref{"ref_2","GGGGGGGGGGGGGGGGGGGGGGGGTAAAAAAAAAAAAAAAAAAAAAAAAG"}
	ref3 := &header_ref{"ref_3",""}
	should_be = append(should_be, ref1,ref2,ref3)

	if len(test_ref) != len(should_be){
		t.Error("Wrong no of refs in test_ref.fa")
	}
	for i :=0; i<len(test_ref); i++ {
		a:=test_ref[i].header
		b:= should_be[i].header
		if a !=b {
			t.Error("Headers dont't match")
		}
		if test_ref[i].seq != should_be[i].seq{
			fmt.Println(test_ref[i].seq)
			fmt.Println(should_be[i].seq)
			t.Error("Seqs dont't match")
		}
	}

}


