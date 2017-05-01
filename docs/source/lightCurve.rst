`LightCurve`
============

The main object representing light curves in `SNData` is the `LightCurve` class. Why might one use the `LightCurve` class to interact with light Curves rather than starting afresh? The idea behind this class (even though the implementation is partial at this point) is that it will:

1. **IO:** Include some IO functionality for reading in data from certain formats that are commonly used. It will also be able to write to certain formats.
2. **Uniformity:** When we read in data from different formats into this, things behave in a uniform way. Thus, apiece of analysis code (for example, fitting a light curve to a model) being written by someone does not need to read such formats and can target this uniform light curve class.
3. **Functionality:**  It provides basic functionality associated with such light curves. For example, it supports `coaddition`, calculating certain simple properties of light curves which can then be used to summarize the light curves.
