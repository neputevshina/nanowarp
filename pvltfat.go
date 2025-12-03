package nanowarp

import (
	"fmt"
	"math"
)

type PVState struct {
	rtpghiState   *RTPGHIState
	stretch       float64
	inPos         int
	outPos        int
	inInOutOffset float64
	outInInOffset float64
	aana          int
	asyn          int
}

func printPos(state *PVState) {
	fmt.Printf("in_pos: %d, out_pos: %d, in_out_offset: %.3f, out_in_offset: %.3f, stretch: %.3f\n",
		state.inPos, state.outPos, state.inInOutOffset, state.outInInOffset, state.stretch)
}

func pvNextInLen(state *PVState, Lout int) int {
	stretch := state.stretch
	inPosEnd := int(math.Round(float64(Lout)/stretch + state.outInInOffset))
	return inPosEnd
}

func pvNextOutLen(state *PVState, Lin int) int {
	stretch := state.stretch
	outPosEnd := int(math.Round(float64(Lin)*stretch + state.inInOutOffset))
	return outPosEnd
}

func pvAdvanceBy(state *PVState, Lin int, Lout int) {
	stretch := state.stretch
	state.inPos += Lin
	state.outPos += Lout
	state.inInOutOffset += float64(Lin) * stretch
	state.inInOutOffset -= float64(Lout)
	state.outInInOffset += float64(Lout) / stretch
	state.outInInOffset -= float64(Lin)
}

func pvSetStretch(state *PVState, stretch float64) int {
	newAana := int(math.Round(float64(state.asyn) / stretch))
	trueStretch := float64(state.asyn) / float64(newAana)

	if trueStretch != state.stretch {
		state.aana = newAana
		state.stretch = trueStretch
	}
	return 0
}

func pvInitStretch(state *PVState, stretch float64) int {
	pvSetStretch(state, stretch)
	rtpghiInitStretch(state.rtpghiState, stretch)
	return 0
}

func rtpghiProcessorCallback(st *PVState, in []complex128, out []complex128) {
	rtpghiExecute(st.rtpghiState, in, st.stretch, out)
}

func pvInit(stretchMax float64, Wmax int, bufLenMax int) (*PVState, error) {
	// asyn := 512
	// M := 8192
	// gl := 4096
	// fifoSize := int(float64(bufLenMax+asyn) * stretchMax)
	// procDelay := gl
	// if gl < fifoSize {
	// 	procDelay = fifoSize
	// }

	// p := &PVState{asyn: asyn}
	// procState, err := RTDGTRProcInitWin(HANN, gl, asyn, M, Wmax, fifoSize, gl)
	// if err != nil {
	// 	return nil, err
	// }
	// p.procState = procState

	// rtpghiState, err := RTPGHIInit(Wmax, asyn, M, 1e-6)
	// if err != nil {
	// 	return nil, err
	// }
	// p.rtpghiState = rtpghiState

	// p.procState.SetCallback(rtpghiProcessorCallback, p)
	// pvSetStretch(p, 1)

	return nil, nil
}

// func pvExecute(p *PVState, in [][]float64, inLen int, chanNo int, stretch float64, outLen int, out [][]float64) int {
// 	pvAdvanceBy(p, inLen, outLen)
// 	pvSetStretch(p, stretch)
// 	// return ExecuteGen(p.procState, in, inLen, chanNo, outLen, out)
// }

// func pvExecuteCompact(p *PVState, in []float64, inLen int, chanNo int, stretch float64, outLen int, out []float64) int {
// 	pvAdvanceBy(p, inLen, outLen)
// 	pvSetStretch(p, stretch)
// 	return p.procState.ExecuteGenCompact(in, inLen, chanNo, outLen, out)
// }

// func pvDone(p **PVState) int {
// 	if p == nil || *p == nil {
// 		return -1
// 	}
// 	pp := *p
// 	if pp.procState != nil {
// 		pp.procState.Done()
// 	}
// 	if pp.rtpghiState != nil {
// 		pp.rtpghiState.Done()
// 	}
// 	*p = nil
// 	return 0
// }
