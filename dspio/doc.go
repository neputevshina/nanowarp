// Experimental package dspio provides primitives and interfaces for basic and granular digital
// signals I/O, including those for granular and overlap-add-based processing.
//
// Please, do not use this package in any contexts not related to Nanowarp or dspio/wavio.
// The API is not standardized and it is not ready for real-time signal processing.
package dspio

// NOTES
// - Do GrainReader/GrainWriter conform to the SignalReader/SignalWriter interface?
//   Probably no, because of ovelaps. Rename those methods to GrainRead/GrainWrite instead.
