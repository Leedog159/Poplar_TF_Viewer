"""Page registry, so pages can navigate to each other with st.switch_page."""
import streamlit as st

from views import about, bulk, orthogroup, single_tf

SINGLE = st.Page(single_tf.render, title="TF Viewer", icon=":material/biotech:",
                 url_path="tf-viewer", default=True)
ORTHO = st.Page(orthogroup.render, title="Orthogroup Viewer", icon=":material/account_tree:",
                url_path="orthogroup")
BULK = st.Page(bulk.render, title="Bulk Lookup & Export", icon=":material/table_view:", url_path="bulk")
ABOUT = st.Page(about.render, title="About", icon=":material/info:", url_path="about")
PAGES = [SINGLE, ORTHO, BULK, ABOUT]
