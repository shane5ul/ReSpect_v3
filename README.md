# ReSpect_v3

An LLM assisted merger of [pyReSpect-freq](https://github.com/shane5ul/pyrespect-freq) and [ReSpect_v2](https://www.mathworks.com/matlabcentral/fileexchange/54322-respect-v2-0).

**Background**

Since about 2016, I basically switched from coding in Matlab/Octave to Python. As I rewrote ReSpect (Matlab) in Python, I made **significant algorithmic improvements** which not only reduced computational time, but also improved quality of the inference. 

These projects are better maintained, and I strongly encourage users of `ReSpect` to use them instead. They are available at:

* `pyrespect-freq`: [$G^{*}(\omega) \rightarrow h(\tau)$](\url{https://github.com/shane5ul/pyReSpect-freq})
* `pyrespect-time`: [$G(t) \rightarrow h(\tau)$](\url{https://github.com/shane5ul/pyReSpect-time})

Consequently, the Matlab code ReSpect-v2 languished over the past decade. With frustration, I noticed that people still used it regularly. It probably has seen more use since it became outdated. Unfortunately, I did not have the bandwidth to maintain ReSpect.

That changed recently. By 2023, AI and LLMs had become quite good at code translation. At the back of my mind, I always thought I might be able to translate my Python projects back to Matlab. This release which used chatGPT, Claude, and Google Gemini to rewrite functions is a result of these actions.

The result is code that is more reliable, and about 10x faster than ReSpect v2. I do not expect to maintain this in the future, but all users of ReSpect_v2.0 are advised to use this version for better performance and accuracy.
