import re

with open('courses/A09.qmd', 'r') as f:
    content = f.read()

# Change eval: false to eval: true and add cache: true and error: true
content = re.sub(r'execute:\n\s*eval:\s*false', 'execute:\n  eval: true\n  cache: true\n  error: true', content)

# Change ```r to ```{r}
content = re.sub(r'^```r$', '```{r}', content, flags=re.MULTILINE)

# Change ```python to ```{python}
content = re.sub(r'^```python$', '```{python}', content, flags=re.MULTILINE)

# Change ```stan to ```{stan, output.var="stan_model"}
content = re.sub(r'^```stan$', '```{stan, output.var="stan_model"}', content, flags=re.MULTILINE)

with open('courses/A09.qmd', 'w') as f:
    f.write(content)

with open('courses/A09.md', 'w') as f:
    f.write(content)

print("Enabled evaluation and updated code block tags.")
