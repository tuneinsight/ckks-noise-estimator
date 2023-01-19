module github.com/tuneinsight/ckks-bootstrapping-precision

go 1.19

require github.com/tuneinsight/lattigo/v4 v4.1.1-0.20230118152750-950e429680cd

replace github.com/tuneinsight/lattigo/v4 => ../lattigo

require (
	golang.org/x/crypto v0.0.0-20220926161630-eccd6366d1be // indirect
	golang.org/x/sys v0.0.0-20220928140112-f11e5e49a4ec // indirect
)
