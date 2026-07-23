.PHONY: test environment-check

test:
	./scripts/sage-python -m pytest

environment-check:
	./scripts/sage-python -m pytest tests/test_environment.py
