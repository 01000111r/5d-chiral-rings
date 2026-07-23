.PHONY: test environment-check

test:
	sage -python -m pytest

environment-check:
	sage -python -m pytest tests/test_environment.py
