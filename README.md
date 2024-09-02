# ReSpect_v3

An LLM assisted merger of [pyReSpect-freq](https://github.com/shane5ul/pyrespect-freq) and [ReSpect_v2](https://www.mathworks.com/matlabcentral/fileexchange/54322-respect-v2-0).

**Background**

Around 2015-2016, I switched my workflows from Matlab/Octave to Python. As a part of this transition, I reimplemented ReSpect (Matlab) in Python while making **significant algorithmic improvements**. As a result pyReSpect is much faster (~10x) and more reliable (using default parameter settings) than the original ReSpect.

These projects are actively maintained, and I strongly encourage users of `ReSpect` to use them instead. They are available at:

* `pyrespect-freq`: [$G^{*}(\omega) \rightarrow h(\tau)$](\url{https://github.com/shane5ul/pyReSpect-freq})
* `pyrespect-time`: [$G(t) \rightarrow h(\tau)$](\url{https://github.com/shane5ul/pyReSpect-time})

Consequently, the Matlab program [ReSpect-v2]((https://www.mathworks.com/matlabcentral/fileexchange/54322-respect-v2-0)) languished over the past decade. Nevertheless, people continued to use it: ReSpect saw more use after it became outdated than before. I watched with frustration as I did not have the will or bandwidth to maintain that project.

That changed recently. By 2023, AI and LLMs had become quite good at code translation. At the back of my mind, I always thought I might be able to translate my Python projects back to Matlab. This particular release which used chatGPT, Claude, and Google Gemini to rewrite functions is the result of this endeavor.

ReSpect_v3 is faster and more reliable than ReSpect v2. I do not expect to maintain this in the future, but all users of ReSpect_v2.0 are advised to use this version for better performance and accuracy.
