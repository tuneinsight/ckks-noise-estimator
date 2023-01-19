module github.com/tuneinsight/ckks-bootstrapping-precision

go 1.19

require (
	github.com/stretchr/testify v1.8.0
	github.com/tuneinsight/lattigo/v4 v4.1.1-0.20230118152750-950e429680cd
)

replace github.com/tuneinsight/lattigo/v4 => ../lattigo

require (
	github.com/davecgh/go-spew v1.1.1 // indirect
	github.com/pmezard/go-difflib v1.0.0 // indirect
	golang.org/x/crypto v0.0.0-20220926161630-eccd6366d1be // indirect
	golang.org/x/sys v0.0.0-20220928140112-f11e5e49a4ec // indirect
	gopkg.in/yaml.v3 v3.0.1 // indirect
)
