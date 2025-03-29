# CollabDrug

A decentralized drug discovery platform that combines blockchain technology with AI for collaborative drug development.

## Project Structure

```
collabdrug/
├── backend/           # FastAPI backend
├── frontend/         # React frontend
├── smart_contracts/  # Solana smart contracts
├── tests/           # Test files
└── docs/            # Documentation
```

## Setup Instructions

### Prerequisites

- Python 3.9+
- Node.js 16+
- PostgreSQL
- Solana CLI tools

### Backend Setup

1. Create a virtual environment:
```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

2. Install dependencies:
```bash
cd backend
pip install -r requirements.txt
```

3. Set up environment variables:
```bash
cp .env.example .env
# Edit .env with your configuration
```

4. Initialize the database:
```bash
# Create database
createdb collabdrug

# Run migrations (to be implemented)
```

5. Start the development server:
```bash
uvicorn api.main:app --reload
```

### Frontend Setup

1. Install dependencies:
```bash
cd frontend
npm install
```

2. Start the development server:
```bash
npm start
```

## Development

- Backend API documentation: http://localhost:8000/docs
- Frontend development server: http://localhost:3000

## Testing

Run tests:
```bash
pytest
```

## Contributing

1. Fork the repository
2. Create a feature branch
3. Commit your changes
4. Push to the branch
5. Create a Pull Request

## License

MIT License 