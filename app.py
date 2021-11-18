import yfinance as yf
import streamlit as st
import main
import numpy as np
import pandas as pd

st.write("""
# SPACAL Module parameters Calculator


""")

wall = st.sidebar.slider("""Wall of Absorber [mm]""", 0.01, 1.00)
capillar_width = st.sidebar.slider('Capillar Tube Width [mm]', 0.00, 1.00)
air_radius = st.sidebar.slider('Air Tube Half Width [mm]', 0.00, 0.50)
fiber_diameter = st.sidebar.slider('Fiber diameter [mm]', 0.1, 4.0)




absorb = st.sidebar.checkbox('Absorber')
moliere_gr = st.sidebar.checkbox('Moliere Radii [gr/cm^2]')
moliere_cm = st.sidebar.checkbox('Moliere Radii [cm]')
rad_len_gr = st.sidebar.checkbox('Radiation Length [gr/cm^2]')
rad_len_cm = st.sidebar.checkbox('Radiation Length [cm]')
all_df = st.sidebar.checkbox('Show the Table')

# 
if st.sidebar.button('Run Calculations'):

    data = main.SummaryFunction(
                                dataframe           =True, 
                                wall                =wall,  # mm   
                                capillar_width      =capillar_width,  # mm
                                air_radius          =air_radius,   # mm 
                                fiber_diameter      =fiber_diameter,   # mm 
                                )

    if absorb:
        # st.write(['st', 'is <', 3])
        st.write("## `Absorber Dimentions`")
        st.text(f'X: {data["X"][0]:.0f} {data["X"][1]}')
        st.text(f'Y: {data["Y"][0]:.0f} {data["Y"][1]}')
        st.text(f'Z: {data["Z"][0]:.0f} {data["Z"][1]}')
    
    if moliere_gr:
        capillar_RM_gr = data["Capillar Moliere Radius [gr/cm^2]"]
        absorber_RM_gr = data["Absorber Moliere Radius [gr/cm^2]"]
        module_RM_gr = data["Module Moliere Radius [gr/cm^2]"]

        st.write("## `Moliere Radii [gr/cm^2]`")
        # st.write("#### RM ")
        st.text(f'Capillar RM: {capillar_RM_gr[0]:.0f} {capillar_RM_gr[1]}')
        st.text(f'Absorber RM: {absorber_RM_gr[0]:.0f} {absorber_RM_gr[1]}')
        st.text(f'Module RM: {module_RM_gr[0]:.0f} {module_RM_gr[1]}')



    if moliere_cm:
        
        capillar_RM_cm = data["Capillar Moliere Radius [cm]"]
        absorber_RM_cm = data["Absorber Moliere Radius [cm]"]
        module_RM_cm = data["Module Moliere Radius [cm]"]


        st.write("## `Moliere Radii  [cm]`")
        # st.write("#### RM")
        st.text(f'Capillar RM: {capillar_RM_cm[0]:.0f} {capillar_RM_cm[1]}')
        st.text(f'Absorber RM: {absorber_RM_cm[0]:.0f} {absorber_RM_cm[1]}')
        st.text(f'Module RM: {module_RM_cm[0]:.0f} {module_RM_cm[1]}')


    
    if rad_len_gr:
        capillar_RL_gr = data["Capillar Rad Length [gr/cm^2]"]
        absorber_RL_gr = data["Absorber Rad Length [gr/cm^2]"]
        module_RL_gr = data["Module Rad Length [gr/cm^2]"]

        st.write("## `Radiation Lengths [gr/cm^2]`")
        st.text(f'Capillar RM: {capillar_RL_gr[0]:.0f} {capillar_RL_gr[1]}')
        st.text(f'Absorber RM: {absorber_RL_gr[0]:.0f} {absorber_RL_gr[1]}')
        st.text(f'Module RM: {module_RL_gr[0]:.0f} {module_RL_gr[1]}')


    if rad_len_cm:
        capillar_RL_cm = data["Capillar Rad Length [cm]"]
        absorber_RL_cm = data["Absorber Rad Length [cm]"]
        module_RL_cm = data["Module Rad Length [cm]"]

        st.write("## `Radiation Lengths [cm]`")
        st.text(f'Capillar RM: {capillar_RL_cm[0]:.0f} {capillar_RL_cm[1]}')
        st.text(f'Absorber RM: {absorber_RL_cm[0]:.0f} {absorber_RL_cm[1]}')
        st.text(f'Module RM: {module_RL_cm[0]:.0f} {module_RL_cm[1]}')

    if all_df:    
        st.dataframe(data.T)




 
