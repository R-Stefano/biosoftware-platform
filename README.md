# Biosoftware Platform (EDEn)

The Electroceutical Design Environment (EDEn) - A comprehensive platform for bioelectricity research and therapeutic design.

## 📋 Credits

This platform was developed under the guidance of **Michael Levin, Ph.D.**

The EDEn platform represents a collaborative effort from the Levin Lab at Tufts University, pioneering research in developmental and regenerative biology through the lens of bioelectricity.

For more information:
- [Levin Lab Website](https://as.tufts.edu/biology/levin-lab)

## 🌟 Overview

EDEn is a pioneering bioinformatics platform that bridges the gap between bioelectricity research and therapeutic applications. The emerging field of bioelectricity has revealed crucial roles for ion channels beyond the nervous system, opening new possibilities in regenerative medicine. This platform serves as an essential tool for researchers and clinicians developing biomedical interventions for various conditions including:
- Birth defects
- Cancer
- Traumatic injury
- Regenerative medicine

### Key Features
- Comprehensive database of ion channels and ion pumps
- Extensive library of chemical modulators and their properties
- Expression data for ion channels across 100+ tissue types
- Interactive visualization tools for bioelectric state analysis
- Therapeutic strategy design interface

### Platform Capabilities
- Identify chemical entities that can modify cellular electrical properties
- Analyze tissue-specific ion channel expression patterns
- Design targeted electroceutical interventions
- Map relationships between ion channels and their modulators

## 🏗️ Architecture

### Frontend (`/frontend`)
- Web application built with Angular
- Select Health vs Disease Tissue, Review Ion Channels expression data, and identify modulators
- Deployment details in app.yaml and cloudbuild.yaml

### Backend (`/backend`)
- Python-based REST API
- Containerized with Docker
- Cloud-ready with GCP integration

### Database (`/db`)
- Scripts to build the database for tissues, expression data, and modulators

## 🚀 Getting Started

### Setup Instructions

1. **Database Setup**
```bash
cd db
pip install -r requirements.txt
Look at README.md in db directory for more details
```

2. **Backend Setup**
```bash
cd backend
pip install -r requirements.txt
python main.py
```

3. **Frontend Setup**
```bash
cd frontend
npm install
node app.js
```
