# Known Bugs and How to Solve Them

## Semaphore stuck when using the large_sample option

``
parallel: Warning: Semaphore stuck for 30 seconds. Consider using --semaphoretimeout.
```

Remove everything in your semaphore directory `~/.parallel/semaphore/`.
