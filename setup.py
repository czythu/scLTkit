# Copyright (C) 2024 Zeyu Chen, G-Lab, Tsinghua University

from setuptools import setup, find_packages
import os

os.chdir(os.path.dirname(os.path.abspath(__file__)))
print(os.getcwd())

with open("README.md", "r") as fh: 
	description = fh.read() 


setup( 
	name="scLTkit",
	version="0.0.3",
	author="Zeyu Chen", 
	author_email="chenzy22@mails.tsinghua.edu.cn", 
    packages=find_packages(),
	description="scLT-kit is a toolkit for analyzing single-cell lineage-tracing (LT-scSeq) data.",
	long_description=description, 
	long_description_content_type="text/markdown", 
	url="https://github.com/czythu/scLTkit",
	license='MIT', 
	python_requires='>=3.7', 
	install_requires=[
		"numpy",
    	"pandas",
    	"scipy",
    	"scikit-learn",
    	"seaborn",
    	"matplotlib",
    	"scanpy",
    	"leidenalg"
	]
) 
