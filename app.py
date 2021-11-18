import streamlit as st
import main
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

st.write("""
# SPACAL Module parameters Calculator


""")

wall = st.sidebar.slider("""Wall of Absorber [mm]""", 0.01, 1.00)
capillar_width = st.sidebar.slider('Capillar Tube Width [mm]', 0.00, 1.00)
air_width = st.sidebar.slider('Air Tube Width [mm]', 0.00, 0.50)
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
                                air_width           =air_width,   # mm 
                                fiber_diameter      =fiber_diameter,   # mm 
                                )

    if absorb:

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
        st.text(f'Capillar RM: {capillar_RM_gr[0]:.4f} {capillar_RM_gr[1]}')
        st.text(f'Absorber RM: {absorber_RM_gr[0]:.4f} {absorber_RM_gr[1]}')
        st.text(f'Module RM: {module_RM_gr[0]:.4f} {module_RM_gr[1]}')



    if moliere_cm:
        
        capillar_RM_cm = data["Capillar Moliere Radius [cm]"]
        absorber_RM_cm = data["Absorber Moliere Radius [cm]"]
        module_RM_cm = data["Module Moliere Radius [cm]"]


        st.write("## `Moliere Radii  [cm]`")
        # st.write("#### RM")
        st.text(f'Capillar RM: {capillar_RM_cm[0]:.4f} {capillar_RM_cm[1]}')
        st.text(f'Absorber RM: {absorber_RM_cm[0]:.4f} {absorber_RM_cm[1]}')
        st.text(f'Module RM: {module_RM_cm[0]:.4f} {module_RM_cm[1]}')


    
    if rad_len_gr:
        capillar_RL_gr = data["Capillar Rad Length [gr/cm^2]"]
        absorber_RL_gr = data["Absorber Rad Length [gr/cm^2]"]
        module_RL_gr = data["Module Rad Length [gr/cm^2]"]

        st.write("## `Radiation Lengths [gr/cm^2]`")
        st.text(f'Capillar RL: {capillar_RL_gr[0]:.4f} {capillar_RL_gr[1]}')
        st.text(f'Absorber RL: {absorber_RL_gr[0]:.4f} {absorber_RL_gr[1]}')
        st.text(f'Module RL: {module_RL_gr[0]:.4f} {module_RL_gr[1]}')


    if rad_len_cm:
        capillar_RL_cm = data["Capillar Rad Length [cm]"]
        absorber_RL_cm = data["Absorber Rad Length [cm]"]
        module_RL_cm = data["Module Rad Length [cm]"]

        st.write("## `Radiation Lengths [cm]`")
        st.text(f'Capillar RL: {capillar_RL_cm[0]:.4f} {capillar_RL_cm[1]}')
        st.text(f'Absorber RL: {absorber_RL_cm[0]:.4f} {absorber_RL_cm[1]}')
        st.text(f'Module RL: {module_RL_cm[0]:.4f} {module_RL_cm[1]}')

    if all_df:    
        st.dataframe(data.T)

draw_2d = st.sidebar.button('Draw 2D Image')


if draw_2d:

    draw_data = main.SummaryFunction(
                                    dataframe           =True, 
                                    wall                =wall,  # mm   
                                    capillar_width      =capillar_width,  # mm
                                    air_width           =air_width,   # mm 
                                    fiber_diameter      =fiber_diameter,   # mm 
                                    )
    
    
    fiber_diameter = fiber_diameter
    # air_width = air_radius * 2
    capillar_width = capillar_width

    hole_diameter = fiber_diameter + air_width + capillar_width

    wall = wall

    fiber_radius = fiber_diameter/2
    air_radius = fiber_radius + air_width/2
    capillar_radius = air_radius + capillar_width/2

    num_fibers = int(draw_data["Fibers per Axis"][0])
    pitch = hole_diameter + wall


    plt.figure(figsize=(10, 10), dpi=500)
    rectangle = plt.Rectangle((0, 0), 121, 121, fc='blue', alpha=0.5)

    plt.gca().add_patch(rectangle)

    zero_coor = (wall/2) + (hole_diameter/2)

    move_pitch = (draw_data["X"][0])/num_fibers

    f"Number of Fibers per Axis: {num_fibers}"
    f"Number of Fibers per Module: {num_fibers**2}"
    x_coor = zero_coor
    y_coor = zero_coor

    for x in range(num_fibers):
        
            for y in range(num_fibers):

                capillar_tube = plt.Circle((x_coor, y_coor), radius = capillar_radius, fc='black')

                air_tube = plt.Circle((x_coor, y_coor), radius = air_radius, fc='blue', alpha=0.5)

                circle = plt.Circle((x_coor, y_coor), radius=fiber_radius, fc='green')

                plt.gca().add_patch(capillar_tube)
                plt.gca().add_patch(air_tube)
                plt.gca().add_patch(circle)

                x_coor += move_pitch
                
                if y == num_fibers-1:
                    
                    x_coor = zero_coor
            
                    y_coor += move_pitch
        

    plt.axis('scaled')
    st.pyplot(plt)




 
