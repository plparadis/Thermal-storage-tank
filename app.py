# Author PLParadis
import streamlit as st
from matplotlib import pyplot as plt


TextColor = "#404041"
plt.rcParams['text.color'] = TextColor
plt.rcParams['axes.labelcolor'] = TextColor
plt.rcParams['xtick.color'] = TextColor
plt.rcParams['ytick.color'] = TextColor
fontname = 'Helvetica LT Std'

st.set_page_config(page_title="",
                    page_icon="static/favicon.ico",
                    layout="centered",
                    initial_sidebar_state="auto")
