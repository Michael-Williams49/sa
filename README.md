# Sequence Alignment Web App (sa)

This web app provides sequence alignment using the Needleman-Wunsch Algorithm (NWA) and Smith-Waterman Algorithm (SWA).

## Installation

1. **Clone the repository:**

```bash
git clone https://github.com/Michael-Williams49/sa
cd sa
```

2. **Install dependencies:**

```bash
pip install -r requirements.txt
```

## Usage

1. Run the app:

```bash
gunicorn app:app
```

2. Open in browser:

Navigate to `localhost:8000` (or the corresponding port if it's different).