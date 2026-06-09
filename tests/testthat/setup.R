# Route any plotting done during tests to a null device. Several plotting
# functions auto-render via `print(p)` on their verbose (`quiet = FALSE`)
# path; without an open device R falls back to writing a stray Rplots.pdf
# into the working directory. Opening pdf(NULL) absorbs those renders.
pdf(NULL)
withr::defer(grDevices::dev.off(), teardown_env())
