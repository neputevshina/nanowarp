package oscope

func Dump(err error, topath string) error {
	if err != nil {
		return err
	}
	if !Enable {
		return nil
	}
	for _, w := range watches {
		switch w.elem.(type) {
		case float64:
			err = dumpWaveform[float64](nil, w, topath)
		case float32:
			err = dumpWaveform[float32](nil, w, topath)
		case int:
			err = dumpWaveform[int](nil, w, topath)
		case uint8:
			err = dumpWaveform[uint8](nil, w, topath)

		case []float64:
			err = dumpTexture[float64](nil, w, topath)
		case []float32:
			err = dumpTexture[float32](nil, w, topath)
		case []int:
			err = dumpTexture[int](nil, w, topath)
		case []uint8:
			err = dumpTexture[uint8](nil, w, topath)
		}
		if err != nil {
			return err
		}
	}
	for _, d := range dists {
		_ = dumpDist(nil, d, topath)
	}
	return nil
}
