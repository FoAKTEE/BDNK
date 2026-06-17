**subagents deployment** 
only allowed type: dynamical workflow claude 4.8 ultracode; claude 4.8 max; all others are NOT ALLOWED

**figure validation**
  - VLM pre-inspection: for each target figure, whose numerical results we want to reproduce, use VLM to inspect the figure before generating plotting script
  - VLM post-dual-inspection: after figure generation, use VLM to inspect both figure for comparisons

**figure generation**
* each time you run a new test, if application, generate a figure
* for the same test which the figure is generated repetitively, git commit the old figure and replace it directly, over-generating too many figures for the same task is prohibited, we can always check later in the commit

**numerical relativity**
  - comparison with original data need to proceed in log scale
    - typical incorrect example: `approximately zero'
    - for special cases, use figure validation, joint DAG, original paper description, to identify the physical meaning of the values in the figure. There is a clear distinction between major component and noise. 
* code clarity: the athenak CCM should be a reusable module
    - files in one subfolder: plain-layout is not a good idea
    - functions in one file: plain-layout is not a good idea
    - use phys-agentic-loop for finding all possiblities to extract abstraction and shorten the code, make it reusable similar to the rest of athenak


