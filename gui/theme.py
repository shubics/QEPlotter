"""Page configuration and restrained scientific-workspace styling."""
import streamlit as st

_CSS = """
<style>
    :root {
        --qep-bg: #111418;
        --qep-sidebar: #0D1013;
        --qep-surface: #181C21;
        --qep-surface-raised: #1D2228;
        --qep-border: #343B44;
        --qep-border-soft: #282E35;
        --qep-text: #E7EAED;
        --qep-muted: #A8B0BA;
        --qep-accent: #7893AE;
        --qep-accent-strong: #9DB2C7;
        --qep-danger: #C77676;
    }

    html, body, [class*="css"] {
        font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Arial,
                     Helvetica, sans-serif;
        color: var(--qep-text);
    }

    h1, h2, h3, h4, h5 {
        color: var(--qep-text) !important;
        font-weight: 600 !important;
        letter-spacing: -0.012em;
    }

    h1 {
        margin-bottom: .35rem !important;
        padding-bottom: .55rem;
        border-bottom: 1px solid var(--qep-border);
        font-size: 2rem !important;
        line-height: 1.25 !important;
    }

    h2 { font-size: 1.45rem !important; }
    h3 { font-size: 1.15rem !important; }

    p, label, .stCaption {
        line-height: 1.55;
    }

    code, pre, [data-testid="stCode"] {
        font-family: "SFMono-Regular", Consolas, "Liberation Mono", monospace;
    }

    .stApp {
        background: var(--qep-bg);
    }

    section[data-testid="stSidebar"] {
        background: var(--qep-sidebar);
        border-right: 1px solid var(--qep-border);
    }

    section[data-testid="stSidebar"] h1 {
        border-bottom: 0;
        font-size: 1.35rem !important;
        letter-spacing: .01em;
    }

    section[data-testid="stSidebar"] [role="radiogroup"] {
        gap: .25rem;
    }

    section[data-testid="stSidebar"] [role="radiogroup"] label {
        min-height: 42px;
        padding: .45rem .55rem;
        border: 1px solid transparent;
        border-radius: 4px;
    }

    section[data-testid="stSidebar"] [role="radiogroup"] label:hover {
        background: var(--qep-surface);
        border-color: var(--qep-border-soft);
    }

    div[data-testid="stMetric"] {
        background: var(--qep-surface);
        border: 1px solid var(--qep-border);
        border-radius: 4px;
        padding: .75rem .9rem;
        box-shadow: none;
    }

    div[data-testid="stMetric"] label {
        color: var(--qep-muted);
    }

    div[data-testid="stMetricValue"] {
        font-variant-numeric: tabular-nums;
    }

    div[data-testid="stFileUploaderDropzone"] {
        background: var(--qep-surface);
        border: 1px dashed #505965;
        border-radius: 4px;
        box-shadow: none;
    }

    button[data-baseweb="tab"] {
        border-radius: 0;
        padding-left: .9rem;
        padding-right: .9rem;
        min-height: 44px;
        color: var(--qep-muted);
    }

    button[data-baseweb="tab"][aria-selected="true"] {
        color: var(--qep-text);
        border-bottom-color: var(--qep-accent);
    }

    div[data-testid="stSlider"] [role="slider"] {
        box-shadow: 0 0 0 3px rgba(120, 147, 174, .20);
    }

    .stButton>button {
        width: 100%;
        min-height: 42px;
        border-radius: 4px;
        border: 1px solid #46505B;
        background: var(--qep-surface-raised);
        color: var(--qep-text);
        font-weight: 550;
        padding: .55rem .8rem;
        box-shadow: none;
        transition: background-color 140ms ease, border-color 140ms ease;
    }

    .stButton>button:hover {
        transform: none;
        background: #242A31;
        border-color: var(--qep-accent);
        color: #FFFFFF;
    }

    button[kind="primary"] {
        background: #4D6882;
        border: 1px solid #6E89A3;
        color: #FFFFFF;
    }

    button[kind="primary"]:hover {
        background: #5B7791;
        border-color: var(--qep-accent-strong);
    }

    button:focus-visible,
    input:focus-visible,
    [role="radio"]:focus-visible,
    [role="tab"]:focus-visible {
        outline: 2px solid var(--qep-accent-strong) !important;
        outline-offset: 2px;
    }

    .stExpander {
        border: 1px solid var(--qep-border) !important;
        border-radius: 4px !important;
        background: var(--qep-surface);
    }

    .stTextInput input,
    .stNumberInput input,
    div[data-baseweb="select"] > div {
        border-radius: 4px !important;
        background: var(--qep-surface) !important;
    }

    div[data-testid="stDataFrame"] {
        border: 1px solid var(--qep-border);
        border-radius: 3px;
        overflow: hidden;
    }

    div[data-testid="stAlert"] {
        border-radius: 4px;
        box-shadow: none;
    }

    hr {
        border-color: var(--qep-border) !important;
    }

    .sym-result {
        margin: .7rem 0 1rem;
        padding: 1rem 1.15rem;
        border: 1px solid var(--qep-border);
        border-left: 3px solid var(--qep-accent);
        border-radius: 4px;
        background: var(--qep-surface);
    }

    .sym-kicker {
        display: block;
        margin-bottom: .45rem;
        color: var(--qep-accent-strong);
        font-size: .72rem;
        font-weight: 700;
        letter-spacing: .08em;
    }

    .sym-equation {
        color: var(--qep-text);
        font-family: Georgia, "Times New Roman", serif;
        font-size: clamp(1.12rem, 2vw, 1.45rem);
        font-weight: 500;
        line-height: 1.45;
        overflow-wrap: anywhere;
    }

    .sym-caption {
        margin-top: .45rem;
        color: var(--qep-muted);
        font-size: .86rem;
        line-height: 1.5;
    }

    @media (prefers-reduced-motion: reduce) {
        *, *::before, *::after {
            scroll-behavior: auto !important;
            transition-duration: 0.01ms !important;
            animation-duration: 0.01ms !important;
            animation-iteration-count: 1 !important;
        }
    }
</style>
"""


def configure_page():
    """Set Streamlit page config (call once, first Streamlit command)."""
    st.set_page_config(
        page_title="QEPlotter",
        page_icon="Q",
        layout="wide",
        initial_sidebar_state="expanded",
    )


def inject_css():
    """Inject the shared dark-theme CSS."""
    st.markdown(_CSS, unsafe_allow_html=True)
