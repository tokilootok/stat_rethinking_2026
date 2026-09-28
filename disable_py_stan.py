import re

def disable_non_r_eval(filepath):
    with open(filepath, 'r') as f:
        content = f.read()

    # Change ```{python} back to ```python
    content = re.sub(r'^```\{python\}$', '```python', content, flags=re.MULTILINE)
    
    # Change ```{python ...} back to ```python
    content = re.sub(r'^```\{python.*?\}$', '```python', content, flags=re.MULTILINE)
    
    # Change ```{stan, ...} back to ```stan
    content = re.sub(r'^```\{stan.*?\}$', '```stan', content, flags=re.MULTILINE)

    with open(filepath, 'w') as f:
        f.write(content)

disable_non_r_eval('courses/A09.qmd')
disable_non_r_eval('courses/A09.md')

print("Disabled evaluation for python and stan blocks.")
