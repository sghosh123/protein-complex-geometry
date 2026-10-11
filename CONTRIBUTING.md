# Contributing

1. Create a branch for a focused change.
2. Add or update synthetic tests for all numerical behavior.
3. Keep coordinate and angle units explicit in names and documentation.
4. Do not commit research trajectories, private structures, absolute paths, cache files, or workstation-specific outputs.
5. Run `python -m pytest -q` and `python -m compileall -q src tests` before opening a pull request.
6. For alignment changes, test a known rigid transform, reflection handling, missing atoms, and non-collinear selections.

Please report expected behavior and include a minimal non-sensitive example for bugs.
