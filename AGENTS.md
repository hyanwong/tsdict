This AGENTS file contains high-level instructions for AI coding assistants working in this repository.

Collaboration / concurrent edits
--------------------------------

- The user may be editing files while an assistant is working; do **not** revert or overwrite unrelated local changes, unless specifically asked.
- Avoid reformatting or touching files you don’t need for the task.

Project structure (top level)
-----------------------------

Tbd as project is built

Tooling overview
----------------

Tbd as project is implemented


Task timeouts (avoid accidental kills)
-------------------------------------

- Set timeouts conservatively for long-running tasks so work is not killed accidentally when it likely would have completed fine.
- For very long runs, it's OK to start them as a background task, `sleep` for a while, check progress (logs/output/files), then `sleep` again until completion.

Coding style and quality
------------------------

- Prioritize semantic correctness, well-defined invariants, and clear error handling over micro-optimizations or clever tricks.
- Structure code to be modular, with coherent responsibilities, good separation of concerns, and an architecture that can be extended without large-scale rewrites.
- Apply DRY carefully: factor out shared logic when it improves clarity and maintenance, but avoid over-abstraction or generic frameworks that obscure intent.
- Aim for readable, explicit code with simple control flow and minimal “magic”; prefer straightforward implementations to compact but opaque idioms.
- Consider algorithmic time and memory complexity when choosing data structures and approaches, especially for code on critical paths or large problem sizes.
- Keep implementations concise but not cryptic; remove dead or unused code and avoid unnecessary indirection.
- For tests, prefer simple, direct assertions that fail fast; test code should not be defensive, over-general, or attempt to mask or recover from failures.
- Scripts should use tqdm (assume installed) to show progress
- Where there is likely to be a long running process and it can be parallelised, use multiprocessing with queues to distribute work across all available CPU cores.
