from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from .routes import molecules, auth

app = FastAPI(
    title="CollabDrug API",
    description="API for decentralized drug discovery platform",
    version="1.0.0"
)

# Configure CORS
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:3000"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Include routers
app.include_router(auth.router, prefix="/api/auth", tags=["auth"])
app.include_router(molecules.router, prefix="/api/molecules", tags=["molecules"])

@app.get("/")
async def root():
    return {"message": "Welcome to CollabDrug API"}

@app.get("/health")
async def health_check():
    return {"status": "healthy"} 