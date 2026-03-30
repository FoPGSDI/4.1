**Managing Skill**
* From https://wxu26.github.io/writings/claude_for_research.html extract two md files for templates. Use these as guiding principles — use the two md files directly, modify them to fit the tasks below, and keep track of progress.

**Central Task**
* Central target: reproduce the multipole expansion of the cow in `arXiv-2504.00506v1/`. Use `arXiv-2504.00506v1/md_output/` as the primary reference for deep research. The final output should be presented in separate sub-folders; make sure code and results are separated. Run plan mode first, then edit mode.
* First call (plan mode): deploy one agent per area —
  - One agent to converge on the convention for mathematical expressions
  - One agent for mathematical non-step-skipping derivations
  - One agent for numerical implementations
  - One agent for test results and test suite design
* Based on the plan, proceed as follows:
* Second call (editing mode): deploy one agent per area (same four as above), each producing code and documentation.
* In editing mode, deploy 11 agents to repeat the task and cross-check.
* Final agent: check and finalize the results.

**Files and Links**
* `arXiv-2504.00506v1/` contains the source paper and converted markdown
* https://github.com/FoPGSDI/4.1 is the git repo, `dev` branch

**Commit Principles**
* Each sub-agent commits individually.
* One commit should be finalized after each stage.
* The progress of each sub-agent should be documented in md files; the final converged convention should be documented in one md file; stage progress should be documented in one md file following the original guideline. These should live in a `progress/` folder.
* We should have stage-wise research diaries and a final report for human reading.
