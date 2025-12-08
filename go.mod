module github.com/neputevshina/nanowarp

go 1.24.0

toolchain go1.24.4

replace github.com/youpy/go-wav => ./wav

require (
	github.com/youpy/go-wav v0.3.2
	golang.org/x/exp v0.0.0-20251002181428-27f1f14c8bb9
	gonum.org/v1/gonum v0.12.1-0.20230202050800-15d6c0edaa65
)

require (
	github.com/youpy/go-riff v0.1.0 // indirect
	github.com/zaf/g711 v1.4.0 // indirect
)
