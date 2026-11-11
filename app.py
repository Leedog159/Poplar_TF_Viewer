"""Poplar TF Activation Domain Atlas.  Run:  streamlit run app.py"""
import streamlit as st

st.set_page_config(page_title="Poplar TF Activation Domain Atlas", page_icon="🌲", layout="wide")

from views.common import inject_css, render_header  # noqa: E402  (after set_page_config)
from views.nav import PAGES  # noqa: E402

inject_css()
nav = st.navigation(PAGES)
render_header()
nav.run()
