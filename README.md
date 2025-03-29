# CollabDrug

<div align="center">
  <img src="assets/images/logo.png" alt="CollabDrug Logo" width="300">
  <p><em>Collaborative Drug Discovery Platform</em></p>
</div>

A decentralized platform for collaborative drug discovery data sharing and collaboration.

## Project Information

- **Website**: [www.collabdrug.xyz](https://www.collabdrug.xyz)
- **Twitter**: [@CollabDrug](https://x.com/CollabDrug)
- **GitHub**: [collabdrug/collabdrug](https://github.com/collabdrug/collabdrug)

## Features

### MVP Version
1. Molecular Data Management
   - Upload molecular data
   - View molecular structures
   - Calculate molecular properties
   - Delete molecular data

2. User Authentication
   - User registration
   - User login
   - Protected routes
   - JWT authentication

3. Basic UI/UX
   - Responsive design
   - Navigation system
   - Error handling
   - Loading states

4. Backend Services
   - RESTful API
   - Database integration
   - Molecular calculation service
   - Authentication service

## Getting Started

### Prerequisites
- Python 3.8+
- Node.js 14+
- PostgreSQL
- Solana CLI tools

### Installation

1. Clone the repository
```bash
git clone https://github.com/collabdrug/collabdrug.git
cd collabdrug
```

2. Install backend dependencies
```bash
cd backend
pip install -r requirements.txt
```

3. Install frontend dependencies
```bash
cd frontend
npm install
```

4. Set up environment variables
```bash
cp .env.example .env
# Edit .env with your configuration
```

5. Run database migrations
```bash
cd backend
alembic upgrade head
```

6. Start the development servers
```bash
# Backend
cd backend
uvicorn main:app --reload

# Frontend
cd frontend
npm start
```

## Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of conduct and the process for submitting pull requests.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details. 