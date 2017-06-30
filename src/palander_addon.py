# This source file contains the Palander add-on for Blender, to be used in
# conjunction with the simulation engine palander_engine.cpp and the scripts
# palander_remoterun.sh and palander_scriptbase. All files are available at
# <https://github.com/FilmakademieRnD/Palander>
#
# Copyright (c) 2017 Animationsinstitut, Filmakademie Baden-Wuerttemberg
# <http://research.animationsinstitut.de/>
#
# This project has been realized in the scope of the Media Solution Center
# project, funded by the state of Baden-Wuerttemberg, Germany:
# <http://msc-bw.de/en>
#
# This file makes use of the content-creation software Blender, which is free
# software released under the GNU General Public License Version 3 (GPLv3).
# Therefore, this file is also released under the GNU GPLv3; you can redistribute
# it and/or modify it under the terms of the GNU GPLv3.
#
# The most recent release of Blender can be downloaded at
# <https://www.blender.org/>
#
# DISCLAIMER: Palander and Blender are free software.
#
# Palander is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the GNU GPLv3 for more details.
#
# You may have received a copy of the GNU General Public License Version 3
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

bl_info = {
    "name": "Palander",
    "description": "Runs the simulation setup on an external engine.",
    "author": "Ari Harjunmaa",
    "version": (1, 0, 1),
    "blender": (2, 78, 0),
    "location": "View3D",
    "support": "COMMUNITY",
    "category": "Import-Export",
}

import bpy
import os
import re
import csv
import subprocess
import stat
import math
from bpy.types import Operator, Panel
from bpy.props import *

def mainCluster(context):
    print('Preparing simulation parameters...')
    os.chdir(bpy.path.abspath("//"))
    with open("toPalabos.xml","w") as params:
        params.write("<?xml ?>\n<!-- Palander configuration file -->\n\n")
        params.write("<Numerics>\n")
        params.write("    <viscosity> " + bpy.context.scene.fluid + "e-6 </viscosity>\n")
        params.write("    <maxTime> " + str(round(bpy.context.scene.simutime,4)) + " </maxTime>\n")
        params.write("    <resolution> " + str(int(round(bpy.context.scene.resolution))) + " </resolution>\n")
        params.write("    <saveOffset> " + str(int(round(bpy.context.scene.saveOffset))) + " </saveOffset>\n")
        params.write("    <latticeVel> " + str(round(bpy.context.scene.latticeVel,4)) + " </latticeVel>\n")
        params.write("    <disableST> " + str(bpy.context.scene.disableST).lower() + " </disableST>\n")
        params.write("    <continuous> " + str(bpy.context.scene.continuous).lower() + " </continuous>\n")
        params.write("    <extraSmooth> " + str(bpy.context.scene.smooooth).lower() + " </extraSmooth>\n")
        params.write("    <useMask> " + str(bpy.context.scene.useMask).lower() + " </useMask>\n")
        params.write("    <nodes> " + str(int(round(bpy.context.scene.nodes))) + " </nodes>\n")
        params.write("    <runhrs> " + str(int(round(bpy.context.scene.runhrs))) + " </runhrs>\n")
        params.write("    <runmins> " + str(int(round(bpy.context.scene.runmins))) + " </runmins>\n")
        params.write("</Numerics>\n\n<Output>\n")
        if bpy.context.scene.useMask:
            params.write("    <maskFile> " + bpy.path.basename(bpy.data.filepath)[:-5] + str(int(round(bpy.context.scene.resolution))).zfill(4) + ".mask </maskFile>\n")
        params.write("    <basePath> " + bpy.path.abspath("//") + " </basePath>\n")
        params.write("    <toCluster> true </toCluster>\n")
        params.write("    <BOBJpre> " + str(bpy.context.scene.BOBJpre).lower() + " </BOBJpre>\n")
        params.write("    <BOBJfin> " + str(bpy.context.scene.BOBJfin).lower() + " </BOBJfin>\n")
        params.write("    <STLint> " + str(bpy.context.scene.STLint).lower() + " </STLint>\n")
        params.write("    <STLsmint> " + str(bpy.context.scene.STLsmint).lower() + " </STLsmint>\n")
        params.write("    <STLmov> " + str(bpy.context.scene.STLmov).lower() + " </STLmov>\n")
        params.write("    <VTKvf> " + str(bpy.context.scene.VTKvf).lower() + " </VTKvf>\n")
        params.write("    <VTKvel> " + str(bpy.context.scene.VTKvel).lower() + " </VTKvel>\n")
        params.write("    <VTKvort> " + str(bpy.context.scene.VTKvort).lower() + " </VTKvort>\n")
        params.write("</Output>\n\n<Geometry>\n")
        oct = 0
        oldFrame = bpy.context.scene.frame_current
        bpy.context.scene.frame_set(bpy.context.scene.saveOffset)
        for obj in bpy.data.objects:
            for mod in obj.modifiers:
                if mod.type == 'FLUID_SIMULATION':
                    params.write("    <object" + str(oct) + ">\n")
                    params.write("        <type" + str(oct) + "> " + mod.settings.type.lower() + " </type" + str(oct) + ">\n")
                    if mod.settings.type == 'DOMAIN':
                        if bpy.context.scene.simuscale == 0.0:
                            params.write("        <size> " + str(round(mod.settings.simulation_scale,4)) + " </size>\n")
                        else:
                            params.write("        <size> " + str(round(bpy.context.scene.simuscale,4)) + " </size>\n")
                        params.write("        <location> " + ' '.join(str(round(obj.location[i],4)) for i in range(0,3)) + " </location>\n")
                        params.write("        <rotation> " + ' '.join(str(round(obj.rotation_euler[i],4)) for i in range(0,3)) + " </rotation>\n")
                        params.write("        <scale> " + ' '.join(str(round(obj.scale[i],4)) for i in range(0,3)) + " </scale>\n")
                        params.write("        <cache> " + mod.settings.filepath + " </cache>\n")
                    elif mod.settings.type == 'INFLOW':
                        params.write("        <inflowVelocity> " + ' '.join(str(round(mod.settings.inflow_velocity[i],4)) for i in range(0,3)) + " </inflowVelocity>\n")
                        if obj.animation_data:
                            if obj.animation_data.action.fcurves[0].data_path == 'modifiers["Fluidsim"].settings.use':
                                params.write("        <enabledFrame> " + ' '.join(str(int(points.co[0])) for points in obj.animation_data.action.fcurves[0].keyframe_points) + " </enabledFrame>\n")
                                params.write("        <enabledSwitch> " + ' '.join(str(int(points.co[1])) for points in obj.animation_data.action.fcurves[0].keyframe_points) + " </enabledSwitch>\n")
                    elif mod.settings.type == 'FLUID':
                        params.write("        <initialVelocity> " + ' '.join(str(round(mod.settings.initial_velocity[i],4)) for i in range(0,3)) + " </initialVelocity>\n")
                    elif mod.settings.type == 'OBSTACLE':
                        if obj.animation_data:
                            cct = int(obj.animation_data.action.fcurves[0].range()[1])
                            params.write("        <motionPathX> " + ' '.join(str(round(obj.animation_data.action.fcurves[0].evaluate(i),4)) for i in range(int(round(bpy.context.scene.saveOffset)),cct)) + " </motionPathX>\n")
                            params.write("        <motionPathY> " + ' '.join(str(round(obj.animation_data.action.fcurves[1].evaluate(i),4)) for i in range(int(round(bpy.context.scene.saveOffset)),cct)) + " </motionPathY>\n")
                            params.write("        <motionPathZ> " + ' '.join(str(round(obj.animation_data.action.fcurves[2].evaluate(i),4)) for i in range(int(round(bpy.context.scene.saveOffset)),cct)) + " </motionPathZ>\n")
                            params.write("        <rotationPathX> " + ' '.join(str(round(obj.animation_data.action.fcurves[3].evaluate(i),4)) for i in range(int(round(bpy.context.scene.saveOffset)),cct)) + " </rotationPathX>\n")
                            params.write("        <rotationPathY> " + ' '.join(str(round(obj.animation_data.action.fcurves[4].evaluate(i),4)) for i in range(int(round(bpy.context.scene.saveOffset)),cct)) + " </rotationPathY>\n")
                            params.write("        <rotationPathZ> " + ' '.join(str(round(obj.animation_data.action.fcurves[5].evaluate(i),4)) for i in range(int(round(bpy.context.scene.saveOffset)),cct)) + " </rotationPathZ>\n")
                    params.write("        <m" + str(oct) + "> ")
                    for vect in obj.matrix_basis:
                        params.write(' '.join(str(round(vect[i],4)) for i in range(0,4)) + " ")
                    params.write("</m" + str(oct) + ">\n        <v" + str(oct) + "> ")
                    for vert in obj.data.vertices:
                        params.write(' '.join(str(round(vert.co[i],4)) for i in range(0,3)) + " ")
                    params.write("</v" + str(oct) + ">\n    </object" + str(oct) + ">\n")
                    if (mod.settings.type == 'OBSTACLE' or mod.settings.type == 'FLUID'):
                        bpy.ops.object.select_all(action='DESELECT')
                        obj.select = True
                        fileName = 'object' + str(oct) + '.stl'
                        if bpy.app.version < (2, 77, 0):
                            bpy.ops.export_mesh.stl(filepath=fileName)
                        else:
                            bpy.ops.export_mesh.stl(filepath=fileName, use_selection=True)
                    oct = oct + 1
        bpy.context.scene.frame_set(oldFrame)
        params.write("    <objectCount> " + str(oct) + " </objectCount>\n")
        params.write("</Geometry>\n")
    print('Submitting simulation to cluster...')
    os.system("scp -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null object*.stl toPalabos.xml " + bpy.context.scene.username + "@hazelhen.hww.de:palander 2>/dev/null")
    os.system("ssh -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null " + bpy.context.scene.username + "@hazelhen.hww.de palander/palander_remoterun.sh 2>/dev/null")
    print('Cleaning up...')
    purge("./", "object[0-9]{1,}.stl")
    os.remove('toPalabos.xml')
    print('Done! Please log on to the cluster to view progress and access simulation results.')

def mainClusterGet(context):
    print('Retrieving simulation results from Hazel Hen...')
    for obj in bpy.data.objects:
        for mod in obj.modifiers:
            if mod.type == 'FLUID_SIMULATION' and mod.settings.type == 'DOMAIN':
                os.chdir(bpy.path.abspath(mod.settings.filepath))
    address = bpy.context.scene.username + "@hazelhen.hww.de"
    proc = subprocess.Popen(["ssh", "-o", "StrictHostKeyChecking=no", "-o", "UserKnownHostsFile=/dev/null", address, "ws_list", "-Cr", "|", "awk", "'NR", "==", "2", "{print($4)}'", "2>/dev/null"], stdout=subprocess.PIPE)
    (out, err) = proc.communicate()
    os.system("scp -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null " + address + ":" + out.decode('utf-8').rstrip() + "/tmp/*.bobj.gz . 2>/dev/null")
    print('Done!')

def mainLocal(context):
    print('Preparing simulation parameters...')
    os.chdir(bpy.path.abspath(bpy.context.scene.localpath))
    if not os.path.exists('tmp'):
        os.makedirs('tmp')
    with open("toPalabos.xml","w") as params:
        params.write("<?xml ?>\n<!-- Palander configuration file -->\n\n")
        params.write("<Numerics>\n")
        params.write("    <viscosity> " + bpy.context.scene.fluid + "e-6 </viscosity>\n")
        params.write("    <maxTime> " + str(round(bpy.context.scene.simutime,4)) + " </maxTime>\n")
        params.write("    <resolution> " + str(int(round(bpy.context.scene.resolution))) + " </resolution>\n")
        params.write("    <saveOffset> " + str(int(round(bpy.context.scene.saveOffset))) + " </saveOffset>\n")
        params.write("    <latticeVel> " + str(round(bpy.context.scene.latticeVel,4)) + " </latticeVel>\n")
        params.write("    <disableST> " + str(bpy.context.scene.disableST).lower() + " </disableST>\n")
        params.write("    <continuous> " + str(bpy.context.scene.continuous).lower() + " </continuous>\n")
        params.write("    <extraSmooth> " + str(bpy.context.scene.smooooth).lower() + " </extraSmooth>\n")
        params.write("    <useMask> " + str(bpy.context.scene.useMask).lower() + " </useMask>\n")
        params.write("</Numerics>\n\n<Output>\n")
        if bpy.context.scene.useMask:
            params.write("    <maskFile> " + bpy.path.basename(bpy.data.filepath)[:-5] + str(int(round(bpy.context.scene.resolution))).zfill(4) + ".mask </maskFile>\n")
        params.write("    <basePath> " + bpy.path.abspath("//") + " </basePath>\n")
        params.write("    <toCluster> false </toCluster>\n")
        params.write("    <BOBJpre> " + str(bpy.context.scene.BOBJpre).lower() + " </BOBJpre>\n")
        params.write("    <BOBJfin> " + str(bpy.context.scene.BOBJfin).lower() + " </BOBJfin>\n")
        params.write("    <STLint> " + str(bpy.context.scene.STLint).lower() + " </STLint>\n")
        params.write("    <STLsmint> " + str(bpy.context.scene.STLsmint).lower() + " </STLsmint>\n")
        params.write("    <STLmov> " + str(bpy.context.scene.STLmov).lower() + " </STLmov>\n")
        params.write("    <VTKvf> " + str(bpy.context.scene.VTKvf).lower() + " </VTKvf>\n")
        params.write("    <VTKvel> " + str(bpy.context.scene.VTKvel).lower() + " </VTKvel>\n")
        params.write("    <VTKvort> " + str(bpy.context.scene.VTKvort).lower() + " </VTKvort>\n")
        params.write("</Output>\n\n<Geometry>\n")
        oct = 0
        oldFrame = bpy.context.scene.frame_current
        bpy.context.scene.frame_set(bpy.context.scene.saveOffset)
        for obj in bpy.data.objects:
            for mod in obj.modifiers:
                if mod.type == 'FLUID_SIMULATION':
                    params.write("    <object" + str(oct) + ">\n")
                    params.write("        <type" + str(oct) + "> " + mod.settings.type.lower() + " </type" + str(oct) + ">\n")
                    if mod.settings.type == 'DOMAIN':
                        if bpy.context.scene.simuscale == 0.0:
                            params.write("        <size> " + str(round(mod.settings.simulation_scale,4)) + " </size>\n")
                        else:
                            params.write("        <size> " + str(round(bpy.context.scene.simuscale,4)) + " </size>\n")
                        params.write("        <location> " + ' '.join(str(round(obj.location[i],4)) for i in range(0,3)) + " </location>\n")
                        params.write("        <rotation> " + ' '.join(str(round(obj.rotation_euler[i],4)) for i in range(0,3)) + " </rotation>\n")
                        params.write("        <scale> " + ' '.join(str(round(obj.scale[i],4)) for i in range(0,3)) + " </scale>\n")
                        params.write("        <cache> " + mod.settings.filepath + " </cache>\n")
                    elif mod.settings.type == 'INFLOW':
                        params.write("        <inflowVelocity> " + ' '.join(str(round(mod.settings.inflow_velocity[i],4)) for i in range(0,3)) + " </inflowVelocity>\n")
                        if obj.animation_data:
                            if obj.animation_data.action.fcurves[0].data_path == 'modifiers["Fluidsim"].settings.use':
                                params.write("        <enabledFrame> " + ' '.join(str(int(points.co[0])) for points in obj.animation_data.action.fcurves[0].keyframe_points) + " </enabledFrame>\n")
                                params.write("        <enabledSwitch> " + ' '.join(str(int(points.co[1])) for points in obj.animation_data.action.fcurves[0].keyframe_points) + " </enabledSwitch>\n")
                    elif mod.settings.type == 'FLUID':
                        params.write("        <initialVelocity> " + ' '.join(str(round(mod.settings.initial_velocity[i],4)) for i in range(0,3)) + " </initialVelocity>\n")
                    elif mod.settings.type == 'OBSTACLE':
                        if obj.animation_data:
                            cct = int(obj.animation_data.action.fcurves[0].range()[1])
                            params.write("        <motionPathX> " + ' '.join(str(round(obj.animation_data.action.fcurves[0].evaluate(i),4)) for i in range(int(round(bpy.context.scene.saveOffset)),cct)) + " </motionPathX>\n")
                            params.write("        <motionPathY> " + ' '.join(str(round(obj.animation_data.action.fcurves[1].evaluate(i),4)) for i in range(int(round(bpy.context.scene.saveOffset)),cct)) + " </motionPathY>\n")
                            params.write("        <motionPathZ> " + ' '.join(str(round(obj.animation_data.action.fcurves[2].evaluate(i),4)) for i in range(int(round(bpy.context.scene.saveOffset)),cct)) + " </motionPathZ>\n")
                            if len(obj.animation_data.action.fcurves) > 3:
                                params.write("        <rotationPathX> " + ' '.join(str(round(obj.animation_data.action.fcurves[3].evaluate(i),4)) for i in range(int(round(bpy.context.scene.saveOffset)),cct)) + " </rotationPathX>\n")
                                params.write("        <rotationPathY> " + ' '.join(str(round(obj.animation_data.action.fcurves[4].evaluate(i),4)) for i in range(int(round(bpy.context.scene.saveOffset)),cct)) + " </rotationPathY>\n")
                                params.write("        <rotationPathZ> " + ' '.join(str(round(obj.animation_data.action.fcurves[5].evaluate(i),4)) for i in range(int(round(bpy.context.scene.saveOffset)),cct)) + " </rotationPathZ>\n")
                    params.write("        <m" + str(oct) + "> ")
                    for vect in obj.matrix_basis:
                        params.write(' '.join(str(round(vect[i],4)) for i in range(0,4)) + " ")
                    params.write("</m" + str(oct) + ">\n        <v" + str(oct) + "> ")
                    for vert in obj.data.vertices:
                        params.write(' '.join(str(round(vert.co[i],4)) for i in range(0,3)) + " ")
                    params.write("</v" + str(oct) + ">\n    </object" + str(oct) + ">\n")
                    if (mod.settings.type == 'OBSTACLE' or mod.settings.type == 'FLUID'):
                        bpy.ops.object.select_all(action='DESELECT')
                        obj.select = True
                        fileName = 'object' + str(oct) + '.stl'
                        if bpy.app.version < (2, 77, 0):
                            bpy.ops.export_mesh.stl(filepath=fileName)
                        else:
                            bpy.ops.export_mesh.stl(filepath=fileName, use_selection=True)
                    oct = oct + 1
        bpy.context.scene.frame_set(oldFrame)
        params.write("    <objectCount> " + str(oct) + " </objectCount>\n")
        params.write("</Geometry>\n")
        
    print('Initializing simulation (this can take quite a while...)')
    crs = int(round(bpy.context.scene.cores))
    if crs == 1:
        subprocess.Popen("./palander_engine")
    else:
        subprocess.Popen(["mpirun", "-np", str(crs), "./palander_engine"])

def mainPresim(context):
    print('Preparing simulation parameters...')
    os.chdir(bpy.path.abspath(bpy.context.scene.localpath))
    if not os.path.exists('out'):
        os.makedirs('out')
    with open("toSPH.xml","w") as params:
        params.write('<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n')
        params.write('<case application="Palander">\n')
        params.write('    <casedef>\n')
        params.write('        <constantsdef>\n')
        params.write('            <lattice bound="1" fluid="1"/>\n')
        params.write('            <gravity x="0.0" y="0.0" z="-9.81"/>\n')
        params.write('            <cflnumber value="0.2"/>\n')
        params.write('            <hswl auto="true" value="0"/>\n')
        params.write('            <coefsound value="10.0"/>\n')
        params.write('            <coefficient value="0.6"/>\n')
        params.write('            <gamma value="7.0"/>\n')
        params.write('            <rhop0 value="1000.0"/>\n')
        params.write('            <eps value="0.5"/>\n')
        params.write('        </constantsdef>\n')
        params.write('        <mkconfig boundcount="230" fluidcount="10"/>\n')
        params.write('        <geometry>\n')
        for obj in bpy.data.objects:
            for mod in obj.modifiers:
                if mod.type == 'FLUID_SIMULATION':
                    if mod.settings.type == 'DOMAIN':
                        if bpy.context.scene.simuscale == 0.0:
                            scale = mod.settings.simulation_scale
                        else:
                            scale = bpy.context.scene.simuscale
                        norm = scale/max(obj.scale)
                        xmin = round((obj.location[0]-obj.scale[0])*norm,4)
                        xmax = round((obj.location[0]+obj.scale[0])*norm,4)
                        ymin = round((obj.location[1]-obj.scale[1])*norm,4)
                        ymax = round((obj.location[1]+obj.scale[1])*norm,4)
                        zmin = round((obj.location[2]-obj.scale[2])*norm,4)
                        zmax = round((obj.location[2]+obj.scale[2])*norm,4)
                        params.write('            <definition dp="' + str(scale/100.0) + '">\n')
                        params.write('                <pointmin x="' + str(xmin) + '" y="' + str(ymin) + '" z="' + str(zmin) + '"/>\n')
                        params.write('                <pointmax x="' + str(xmax) + '" y="' + str(ymax) + '" z="' + str(zmax) + '"/>\n')
                        params.write('            </definition>\n')
        params.write('            <commands>\n')
        params.write('                <mainlist>\n')
        params.write('                    <setshapemode>actual | dp | fluid | bound | void</setshapemode>\n')
        params.write('                    <setdrawmode mode="full"/>\n')
        oct = 0
        flind = -1
        for obj in bpy.data.objects:
            for mod in obj.modifiers:
                if mod.type == 'FLUID_SIMULATION':
                    if (mod.settings.type == 'FLUID' or mod.settings.type == 'OBSTACLE'):
                        bpy.ops.object.select_all(action='DESELECT')
                        obj.select = True
                        fileName = 'object' + str(oct) + '.stl'
                        if bpy.app.version < (2, 77, 0):
                            bpy.ops.export_mesh.stl(filepath=fileName)
                        else:
                            bpy.ops.export_mesh.stl(filepath=fileName, use_selection=True)
                        if mod.settings.type == 'OBSTACLE':
                            if flind == -1:
                                flind = oct
                            params.write('                    <setmkbound B="0" G="0" R="255" mk="' + str(flind) + '" name="' + obj.name + '"/>\n')
                            params.write('                    <drawfilestl file="' + fileName + '">\n')
                            params.write('                        <drawscale x="' + str(0.5*norm) + '" y="' + str(0.5*norm) + '" z="' + str(0.5*norm) + '"/>\n')
                            params.write('                    </drawfilestl>\n')
                        else:
                            params.write('                    <setmkfluid B="255" G="0" R="0" mk="' + str(oct) + '" name="' + obj.name + '"/>\n')
                            params.write('                    <drawfilestl file="' + fileName + '">\n')
                            params.write('                        <drawscale x="' + str(0.5*norm) + '" y="' + str(0.5*norm) + '" z="' + str(0.5*norm) + '"/>\n')
                            params.write('                    </drawfilestl>\n')
                            params.write('                    <fillpoint x="' + str(round(0.5*norm*obj.location[0],4)) + '" y="' + str(round(0.5*norm*obj.location[1],4)) + '" z="' + str(round(0.5*norm*obj.location[2],4)) + '" mkfluid="' + str(oct) + '">\n')
                            params.write('                        <modefill>void</modefill>\n')
                            params.write('                    </fillpoint>\n')
                        params.write('                    <shapeout file="' + obj.name + '"/>\n')
                    oct = oct + 1
        for obj in bpy.data.objects:
            for mod in obj.modifiers:
                if mod.type == 'FLUID_SIMULATION':
                    if mod.settings.type == 'DOMAIN':
                        params.write('                    <setmkbound mk="' + str(oct) + '"/>\n')
                        params.write('                    <drawbox>\n')
                        params.write('                        <boxfill>bottom | left | right | front | back</boxfill>\n')
                        params.write('                        <point x="' + str(round((obj.location[0]-obj.scale[0])*0.5*norm,4)) + '" y="' + str(round((obj.location[1]-obj.scale[1])*0.5*norm,4)) + '" z="' + str(round((obj.location[2]-obj.scale[2])*0.5*norm,4)) +'"/>\n')
                        params.write('                        <size x="' + str(round(obj.scale[0]*norm,4)) + '" y="' + str(round(obj.scale[1]*norm,4)) + '" z="' + str(round(obj.scale[2]*norm,4)) +'"/>\n')
                        params.write('                    </drawbox>\n')
        params.write('                </mainlist>\n')
        params.write('            </commands>\n')
        params.write('        </geometry>\n')
        params.write('        <floatings>\n')
        params.write('            <floating mkbound="' + str(flind) + '" relativeweight="0.8"/>\n') #
        params.write('        </floatings>\n')
        params.write('    </casedef>\n')
        params.write('    <execution>\n')
        params.write('        <parameters>\n')
        params.write('            <parameter key="StepAlgorithm" value="1"/>\n')
        params.write('            <parameter key="VerletSteps" value="40"/>\n')
        params.write('            <parameter key="Kernel" value="1"/>\n')
        params.write('            <parameter key="ViscoTreatment" value="1"/>\n')
        params.write('            <parameter key="Visco" value="0.1"/>\n')
        params.write('            <parameter key="ShepardSteps" value="0"/>\n')
        params.write('            <parameter key="DeltaSPH" value="0"/>\n')
        params.write('            <parameter key="DtIni" value="1.0e-4.0"/>\n')
        params.write('            <parameter key="DtMin" value="1.0e-6.0"/>\n')
        params.write('            <parameter key="TimeMax" value="' + str(int(round(bpy.context.scene.simutime))) + '"/>\n')
        params.write('            <parameter key="TimeOut" value="0.04"/>\n') # 1/FPS
        params.write('            <parameter key="IncZ" value="0.01"/>\n')
        params.write('            <parameter key="PartsOutMax" value="1"/>\n')
        params.write('            <parameter key="RhopOutMax" value="1300.0"/>\n')
        params.write('            <parameter key="RhopOutMin" value="700.0"/>\n')
        params.write('        </parameters>\n')
        params.write('    </execution>\n')
        params.write('</case>\n')
    print('Starting floating pre-simulation...')
    with open("presim.sh","w") as script:
        sphdir = bpy.path.abspath(bpy.context.scene.SPHpath)
        script.write('#!/bin/bash\n\n')
        script.write('gencase="' + sphdir + '/EXECS/GenCase4_linux64"\n')
        script.write('dualsphysics="' + sphdir + '/EXECS/DualSPHysics4CPU_linux64"\n')
        script.write('floatinginfo="' + sphdir + '/EXECS/FloatingInfo4_linux64"\n')
        script.write('path_so="' + sphdir + '/EXECS"\n')
        script.write('export LD_LIBRARY_PATH=$path_so\n')
        script.write('$gencase toSPH out/presim -save:all\n')
        script.write('$dualsphysics out/presim out -svres -cpu\n')
        script.write('$floatinginfo -filexml out/presim.xml -savemotion -savedata out/motionPath\n')
    st = os.stat("presim.sh")
    os.chmod("presim.sh", st.st_mode | stat.S_IEXEC)
    subprocess.Popen(["./presim.sh"])
    print('Done! Please wait for the pre-simulation to complete before accessing the results.')

def mainRetrieve(context):
    print('Retrieving float path results...')
    os.chdir(bpy.path.abspath(bpy.context.scene.localpath))
    if not os.path.exists('out'):
        os.makedirs('out')
    found = False
    for f in os.listdir("out"):
        if re.search("motionPath_mk[0-9]{1,}.csv", f):
            found = True
            with open("out/" + f, "r") as motionfile:
                motion = csv.reader(motionfile, delimiter=';')
                next(motion, None)
                data = list(motion)
    if not found:
        print('Floating data not available!')
    else:
        for obj in bpy.data.objects:
            for mod in obj.modifiers:
                if (mod.type == 'FLUID_SIMULATION' and mod.settings.type == 'DOMAIN'):
                    if bpy.context.scene.simuscale == 0.0:
                        scale = mod.settings.simulation_scale
                    else:
                        scale = bpy.context.scene.simuscale
                    norm = 0.5*scale/max(obj.scale)
                if (mod.type == 'FLUID_SIMULATION' and mod.settings.type == 'OBSTACLE'):
                    radtodeg = math.pi/180.0
                    for i in range (len(data)):
                        obj.location = [float(data[i][8])/norm,float(data[i][9])/norm,float(data[i][10])/norm]
                        obj.rotation_euler = [float(data[i][14])*radtodeg,float(data[i][16])*radtodeg,float(data[i][15])*radtodeg]
                        toframe = int(data[i][0])
                        obj.keyframe_insert(data_path="location",frame=toframe)
                        obj.keyframe_insert(data_path="rotation_euler",frame=toframe)
                    for fc in obj.animation_data.action.fcurves:
                        fc.update()
        print('Float path imported!')

def initProps():
    bpy.types.Scene.fluid = bpy.props.EnumProperty(
        items = [('4',"Blood",'Blood (4)'),
            ('0.6',"Gasoline",'Gasoline (0.6)'),
            ('2000',"Honey",'Honey (2000)'),
            ('50',"Oil",'Oil (50)'),
            ('1',"Water",'Water (1)')],
        name = "Fluid",
        default = '1'
        )

    bpy.types.Scene.simutime = bpy.props.FloatProperty(
        name = "Time",
        default = 5,
        min = 0.1,
        max = 100,
        step = 10,
        precision = 1,
        description = "Simulation length in seconds"
        )

    bpy.types.Scene.simuscale = bpy.props.FloatProperty(
        name = "Size",
        default = 0,
        min = 0,
        max = 1000,
        step = 10,
        precision = 3,
        description = "Real world size in meters (0 to use Blender value)"
        )

    bpy.types.Scene.resolution = bpy.props.FloatProperty(
        name = "Resolution",
        default = 200,
        min = 50,
        max = 2000,
        step = 1000,
        precision = 0,
        description = "Mesh resolution"
        )

    bpy.types.Scene.saveOffset = bpy.props.FloatProperty(
        name = "Frame offset",
        default = 0,
        min = -10000,
        max = 10000,
        step = 100,
        precision = 0,
        description = "Frame at which to start simulation"
        )

    bpy.types.Scene.latticeVel = bpy.props.FloatProperty(
        name = "Lattice velocity",
        default = 0.025,
        min = 0.001,
        max = 0.1,
        step = 0.1,
        precision = 3,
        description = "Lattice velocity of system"
        )

    bpy.types.Scene.disableST = bpy.props.BoolProperty(
        name = "Disable surface tension",
        default = False,
        description = "Disable surface tension algorithm for extra splashiness"
        )

    bpy.types.Scene.continuous = bpy.props.BoolProperty(
        name = "Continuous simulation",
        default = False,
        description = "Simulate continuously using saved checkpoints"
        )

    bpy.types.Scene.smooooth = bpy.props.BoolProperty(
        name = "Extra-smooth surface",
        default = False,
        description = "Produce extra-smooth surface at the cost of detail"
        )

    bpy.types.Scene.useMask = bpy.props.BoolProperty(
        name = "Use Mask",
        default = False,
        description = "Load/save boolean mask"
        )

    bpy.types.Scene.SPHpath = bpy.props.StringProperty(
        name = "Path",
        subtype = "DIR_PATH",
        description = "Path to DualSPHysics root directory"
        )

    bpy.types.Scene.BOBJpre = bpy.props.BoolProperty(
        name = "Preview Surface",
        default = True,
        description = "Output preview-quality bobj.gz surface"
        )

    bpy.types.Scene.BOBJfin = bpy.props.BoolProperty(
        name = "Final Surface",
        default = True,
        description = "Output final-quality bobj.gz surface"
        )

    bpy.types.Scene.STLint = bpy.props.BoolProperty(
        name = "Preview Interface",
        default = True,
        description = "Output preview-quality STL interface"
        )

    bpy.types.Scene.STLsmint = bpy.props.BoolProperty(
        name = "Final Interface",
        default = True,
        description = "Output final-quality, smoothened STL interface"
        )

    bpy.types.Scene.STLmov = bpy.props.BoolProperty(
        name = "Moving Obstacle",
        default = True,
        description = "Output moving obstacle STL files"
        )

    bpy.types.Scene.VTKvf = bpy.props.BoolProperty(
        name = "Volume Fraction",
        default = False,
        description = "Include volume fractions in VTK output"
        )

    bpy.types.Scene.VTKvel = bpy.props.BoolProperty(
        name = "Velocity",
        default = False,
        description = "Include velocity field in VTK output"
        )

    bpy.types.Scene.VTKvort = bpy.props.BoolProperty(
        name = "Vorticity",
        default = False,
        description = "Include vorticity field in VTK output"
        )

    bpy.types.Scene.localpath = bpy.props.StringProperty(
        name = "Path",
        subtype = "DIR_PATH",
        description = "Path to local simulation directory"
        )

    bpy.types.Scene.cores = bpy.props.FloatProperty(
        name = "Cores",
        default = 1,
        min = 1,
        max = 32,
        step = 100,
        precision = 0,
        description = "Number of cores to use on the local machine"
        )

    bpy.types.Scene.username = bpy.props.StringProperty(
        name = "User",
        description = "User account name on Hazel Hen"
        )

    bpy.types.Scene.runhrs = bpy.props.FloatProperty(
        name = "Hr",
        default = 0,
        min = 0,
        max = 23,
        step = 100,
        precision = 0,
        description = "Wall time (hours)"
        )

    bpy.types.Scene.runmins = bpy.props.FloatProperty(
        name = "Min",
        default = 25,
        min = 0,
        max = 60,
        step = 100,
        precision = 0,
        description = "Wall time (minutes)"
        )

    bpy.types.Scene.nodes = bpy.props.FloatProperty(
        name = "Nodes",
        default = 10,
        min = 1,
        max = 7712,
        step = 100,
        precision = 0,
        description = "Number of nodes to use on the cluster"
        )

def removeProps():
    del bpy.types.Scene.fluid
    del bpy.types.Scene.simutime
    del bpy.types.Scene.simuscale
    del bpy.types.Scene.resolution
    del bpy.types.Scene.saveOffset
    del bpy.types.Scene.latticeVel
    del bpy.types.Scene.disableST
    del bpy.types.Scene.continuous
    del bpy.types.Scene.smooooth
    del bpy.types.Scene.useMask
    del bpy.types.Scene.SPHpath
    del bpy.types.Scene.BOBJpre
    del bpy.types.Scene.BOBJfin
    del bpy.types.Scene.STLint
    del bpy.types.Scene.STLsmint
    del bpy.types.Scene.STLmov
    del bpy.types.Scene.VTKvf
    del bpy.types.Scene.VTKvel
    del bpy.types.Scene.VTKvort
    del bpy.types.Scene.localpath
    del bpy.types.Scene.cores
    del bpy.types.Scene.username
    del bpy.types.Scene.runhrs
    del bpy.types.Scene.runmins
    del bpy.types.Scene.nodes

def purge(dir, pattern):
    for f in os.listdir(dir):
        if re.search(pattern, f):
            os.remove(os.path.join(dir, f))

class Parameters(bpy.types.Panel):
    bl_label = "Simulation Parameters"
    bl_idname = "OBJECT_PT_PARAM"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'TOOLS'
    bl_context = "objectmode"
    bl_category = "Palander"

    def draw(self, context):
        layout = self.layout
        scn = context.scene
        col = layout.column(align=True)
        row = col.row(align=True)
        row.prop(scn, "fluid")
        row = col.row(align=True)
        row.prop(scn, "simutime")
        row = col.row(align=True)
        row.prop(scn, "simuscale")
        row = col.row(align=True)
        row.prop(scn, "resolution")
        row = col.row(align=True)
        row.prop(scn, "saveOffset")
        row = col.row(align=True)
        row.prop(scn, "latticeVel")
        row = col.row(align=True)
        row.label(text="")
        row = col.row(align=True)
        row.prop(scn, "disableST")
        row = col.row(align=True)
        row.prop(scn, "continuous")
        row = col.row(align=True)
        row.prop(scn, "smooooth")

class Floating(bpy.types.Panel):
    bl_label = "DualSPHysics"
    bl_idname = "OBJECT_PT_FLOAT"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'TOOLS'
    bl_context = "objectmode"
    bl_category = "Palander"

    def draw(self, context):
        layout = self.layout
        scn = context.scene
        col = layout.column(align=True)
        row = col.row(align=True)
        row.prop(scn, "SPHpath")
        row = col.row(align=True)
        row.operator("object.pre_simulation", text="Floating Pre-Simulation")
        row = col.row(align=True)
        row.operator("object.retrieve", text="Retrieve Float Path")

class HLRSAccess(bpy.types.Panel):
    bl_label = "Palabos"
    bl_idname = "OBJECT_PT_HLRS"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'TOOLS'
    bl_context = "objectmode"
    bl_category = "Palander"

    def draw(self, context):
        layout = self.layout
        scn = context.scene
        col = layout.column(align=True)
        row = col.row(align=True)
        row.prop(scn, "useMask")
        row = col.row(align=True)
        row.label(text="BOBJ Output:")
        row = col.row(align=True)
        row.prop(scn, "BOBJpre")
        row = col.row(align=True)
        row.prop(scn, "BOBJfin")
        row = col.row(align=True)
        row.label(text="STL Output:")
        row = col.row(align=True)
        row.prop(scn, "STLint")
        row = col.row(align=True)
        row.prop(scn, "STLsmint")
        row = col.row(align=True)
        row.prop(scn, "STLmov")
        row = col.row(align=True)
        row.label(text="VTK Output:")
        row = col.row(align=True)
        row.prop(scn, "VTKvf")
        row = col.row(align=True)
        row.prop(scn, "VTKvel")
        row = col.row(align=True)
        row.prop(scn, "VTKvort")
        row = col.row(align=True)
        row.label(text="Local Simulation:")
        row = col.row(align=True)
        row.prop(scn, "localpath")
        row = col.row(align=True)
        row.prop(scn, "cores")
        row = col.row(align=True)
        row.operator("object.local_send", text="Execute Locally")
        row = col.row(align=True)
        row.label(text="Cluster Simulation:")
        row = col.row(align=True)
        row.prop(scn, "username")
        row = col.row(align=True)
        row.prop(scn, "nodes")
        row = col.row(align=True)
        split = row.split(align=True, percentage=0.5)
        split.prop(scn, "runhrs")
        split.prop(scn, "runmins")
        row = col.row(align=True)
        row.operator("object.hlrs_send", text="Send to Cluster")
        row = col.row(align=True)
        row.operator("object.hlrs_get", text="Retrieve Results")

class ClusterSend(bpy.types.Operator):
    """Connect to HPC Cluster"""
    bl_idname = "object.hlrs_send"
    bl_label = "Send to Cluster"
    bl_options = {'REGISTER'}

    @classmethod
    def poll(cls, context):
        return context.active_object is not None

    def execute(self, context):
        mainCluster(context)
        return {'FINISHED'}

class ClusterGet(bpy.types.Operator):
    """Retrieve Results from HPC Cluster"""
    bl_idname = "object.hlrs_get"
    bl_label = "Retrieve Results"
    bl_options = {'REGISTER'}

    @classmethod
    def poll(cls, context):
        return context.active_object is not None

    def execute(self, context):
        mainClusterGet(context)
        return {'FINISHED'}

class LocalSend(bpy.types.Operator):
    """Execute Locally"""
    bl_idname = "object.local_send"
    bl_label = "Execute Locally"
    bl_options = {'REGISTER'}

    @classmethod
    def poll(cls, context):
        return context.active_object is not None

    def execute(self, context):
        mainLocal(context)
        return {'FINISHED'}

class PreSimulation(bpy.types.Operator):
    """Execute Pre-Simulation"""
    bl_idname = "object.pre_simulation"
    bl_label = "Floating Pre-Simulation"
    bl_options = {'REGISTER'}

    @classmethod
    def poll(cls, context):
        return context.active_object is not None

    def execute(self, context):
        mainPresim(context)
        return {'FINISHED'}

class Retrieve(bpy.types.Operator):
    """Retrieve Float Path"""
    bl_idname = "object.retrieve"
    bl_label = "Retrieve Float Path"
    bl_options = {'REGISTER'}

    @classmethod
    def poll(cls, context):
        return context.active_object is not None

    def execute(self, context):
        mainRetrieve(context)
        return {'FINISHED'}

def register():
    bpy.utils.register_module(__name__)
    initProps()

def unregister():
    bpy.utils.unregister_module(__name__)
    removeProps()

if __name__ == "__main__":
    register()
