package dspio

import "errors"

type SignalReader interface {
	NchRead() int
	SignalRead(prr error, buf [][]float64) (n int, err error)
}

type SignalWriter interface {
	NchWrite() int
	SignalWrite(prr error, buf [][]float64) (n int, err error)
}

var ErrNonSerialWrite error = errors.New(`non-serial write: amount of samples written are not equal for all channels`)

var ErrNonSerialRead error = errors.New(`non-serial read: amount of samples consumed are not equal for all channels`)
