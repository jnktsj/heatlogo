#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2012, 2013 Junko Tsuji

# Convert disorder and secondary structure results to
# Sequence Logo like format.


from corebio.seq import Alphabet

from weblogolib.color import *
from weblogolib.colorscheme import HeatScheme, HeatGroup
from weblogolib.colorscheme import ColorScheme, ColorGroup


# symbols of disorder and secondary structure
structure_alphabet = Alphabet("DCHE", zip('dche', 'DCHE'))


# members of heatmap color schemes
thermography = HeatScheme( [
    HeatGroup( '9', '#FF0000' ),
    HeatGroup( '8', '#FF6600' ),
    HeatGroup( '7', '#FFCC00' ),
    HeatGroup( '6', '#FFFF00' ),
    HeatGroup( '5', '#CCFF33' ),
    HeatGroup( '4', '#00FFCC' ),
    HeatGroup( '3', '#00FFFF' ),
    HeatGroup( '2', '#0066FF' ),
    HeatGroup( '1', '#0000FF' ),
    HeatGroup( '0', '#330099' ) ],
    title = "Thermographic colors" )

red_green = HeatScheme( [
    HeatGroup( '9', '#FF0000' ),
    HeatGroup( '8', '#FF2020' ),
    HeatGroup( '7', '#FF4040' ),
    HeatGroup( '6', '#FF8080' ),
    HeatGroup( '5', '#FFC0C0' ),
    HeatGroup( '4', '#C0FFC0' ),
    HeatGroup( '3', '#80FF80' ),
    HeatGroup( '2', '#40FF40' ),
    HeatGroup( '1', '#20FF20' ),
    HeatGroup( '0', '#00FF00' ) ],
    title = "Gradational colors between red and green" )

red_blue = HeatScheme( [
    HeatGroup( '9', '#FF0000' ),
    HeatGroup( '8', '#FF2020' ),
    HeatGroup( '7', '#FF4040' ),
    HeatGroup( '6', '#FF8080' ),
    HeatGroup( '5', '#FFC0C0' ),
    HeatGroup( '4', '#C0C0FF' ),
    HeatGroup( '3', '#8080FF' ),
    HeatGroup( '2', '#4040FF' ),
    HeatGroup( '1', '#2020FF' ),
    HeatGroup( '0', '#0000FF' ) ],
    title = "Gradational colors between red and blue" )

magenta_cyan = HeatScheme( [
    HeatGroup( '9', '#FF00FF' ),
    HeatGroup( '8', '#FF40FF' ),
    HeatGroup( '7', '#FF60FF' ),
    HeatGroup( '6', '#FFA0FF' ),
    HeatGroup( '5', '#FFE0FF' ),
    HeatGroup( '4', '#E0FFFF' ),
    HeatGroup( '3', '#A0FFFF' ),
    HeatGroup( '2', '#60FFFF' ),
    HeatGroup( '1', '#40FFFF' ),
    HeatGroup( '0', '#00FFFF' ) ],
    title = "Gradational colors between magenta and cyan" )


# heatmap color schemes of disorder and secondary structure
physichem_color_schemes = {"thermography": thermography,
                           "red-green": red_green,
                           "red-blue": red_blue,
                           "magenta-cyan": magenta_cyan}
