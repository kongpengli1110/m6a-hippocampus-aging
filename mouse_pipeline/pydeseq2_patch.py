# joblib cannot use multiprocessing in this sandbox (no /dev/shm) -> neutralise parallel_backend in pydeseq2
import contextlib, pydeseq2.default_inference as _di
_di.parallel_backend=lambda *a,**k: contextlib.nullcontext()
