Small upstream candidates extracted from wojtek-dev, 2026-10-24.

These patches contain no GLM implementation or benchmark archive. They were
prepared against PR80 and tested together in an isolated PR80 source archive.
The three focused suites pass 38 assertions. See inst/WOJTEK-DEV.txt for the
full review and validation details. No pull requests have been submitted.

01-covariance.patch: vcov()/coef() visibility and prediction intervals.
02-compressed-cluster.patch: expanded interaction source for clustered scores.
03-sql-aliases.patch: collision-safe aliases, including utils::tail qualification.

On a new branch from the intended upstream base, use git apply --check first,
then git apply for the selected patch and run the package's full test/check CI.
Each patch adds its own regression test; do not include this patch directory
or the rest of wojtek-dev in the submitted PR. If upstream has moved, adapt
only the surrounding context and revalidate.
