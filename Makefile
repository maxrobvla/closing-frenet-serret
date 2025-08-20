all:
	python make

.PHONY: debug
debug:
	python make --debug

.PHONY: default
default:
	python make --default-generator

.PHONY: default_debug
default_debug:
	python make --default-generator --debug
