def is_running_from_ipython():
    from IPython import get_ipython
    return get_ipython() is not None

if is_running_from_ipython():
    import src.gui.buckinghampigui as gui
