package scram2pkg

import (
	"bytes"
	"fmt"
	"sync"
	"time"
)

//AlignReads aligns reads of  length nt to one or more reference sequences, with exact matches in forward or reverse
//complement accepted.
//A map of ref_header:[srna_seq:[pos,pos,...],...] is returned.
func AlignReads(seq_map map[string]*mean_se, ref_slice []*header_ref, nt int) map[string]map[string][]int {
	wg := &sync.WaitGroup{}
	t1 := time.Now()
	wg.Add(len(ref_slice))

	ref_seq_chan := make(chan *header_ref, len(ref_slice))
	for _, header_ref_pair := range ref_slice {
		ref_seq_chan <- header_ref_pair
	}
	close(ref_seq_chan)

	header_map_chan := make(chan map[string]map[string][]int, 1)
	for a := 0; a < len(ref_slice); a++ {
		go worker_go(seq_map, ref_seq_chan, nt, header_map_chan, wg)
	}
	go func(cs chan map[string]map[string][]int, wg *sync.WaitGroup) {
		wg.Wait()
		close(cs)
	}(header_map_chan, wg)

	final_alignment_map := compile_alignments(header_map_chan)
	t2 := time.Since(t1)
	fmt.Println("Read file set aligned to reference/s: ", t2)
	return final_alignment_map
}

//Align reads to individual reference sequences
func worker_go(seq_map map[string]*mean_se, ref_seqs chan *header_ref, nt int,
	header_map_chan chan map[string]map[string][]int, wg *sync.WaitGroup) {

	ref_seq := <-ref_seqs
	rvs_seq := reverse_complement(ref_seq.seq)
	//srna alignment map --> srna_seq:[pos, pos]
	header_mapped := make(map[string][]int)
	//header_alignment map header:[srna_seq:[pos,pos]]
	mapped := make(map[string]map[string][]int, 1)
	position := 0
	ref_seq_len := len(ref_seq.seq)
	for position <= ref_seq_len-nt {
		//each alignment position for an srna
		fwd_seq := ref_seq.seq[position : position+nt]
		rvs_seq := rvs_seq[position : position+nt]
		if _, ok := seq_map[fwd_seq]; ok {
			header_mapped[fwd_seq] = append(header_mapped[fwd_seq], 1+position)
		}
		if _, ok := seq_map[rvs_seq]; ok {
			header_mapped[rvs_seq] = append(header_mapped[rvs_seq], -1-(ref_seq_len-position-nt))
		}
		position++
	}
	if len(header_mapped) > 0 {
		mapped[ref_seq.header] = header_mapped
	}
	header_map_chan <- mapped
	wg.Done()
}

//reverse_complement a DNA sequence
func reverse_complement(dna string) string {
	complement := map[rune]rune{
		'A': 'T',
		'C': 'G',
		'G': 'C',
		'T': 'A',
		'N': 'N',
	}
	runes := []rune(dna)
	var result bytes.Buffer
	for i := len(runes) - 1; i >= 0; i-- {
		result.WriteRune(complement[runes[i]])
	}
	return result.String()
}

//compile_alignments compiles the alignments
func compile_alignments(header_map_chan chan map[string]map[string][]int) map[string]map[string][]int {
	final_alignment_map := make(map[string]map[string][]int)
	for combined_alignment := range header_map_chan {
		for header, alignment := range combined_alignment {
			final_alignment_map[header] = alignment
		}
	}
	return final_alignment_map
}
