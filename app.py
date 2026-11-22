import streamlit as st


st.set_page_config(page_title="Slabber", page_icon=":material/grid_on:", layout="wide")
page = st.navigation([
    st.Page(
        "app_pages/four_edge.py", title="Four-edge slabs", icon=":material/grid_on:", default=True,
        help="Explore the nine Czerny support cases and their finite-difference benchmarks.",
    ),
    st.Page(
        "app_pages/three_edge.py", title="Three-edge slabs", icon=":material/border_bottom:",
        help="Explore the six Hahn free-edge cases and their Rayleigh–Ritz benchmarks.",
    ),
])
with st.sidebar:
    st.subheader("Slabber", anchor=False)
    st.caption("Elastic plate benchmarks")

page.run()
