from sqlalchemy import create_engine
from functools import lru_cache
from pydantic_settings import BaseSettings

class Settings(BaseSettings):
    DB_URL: str = "sqlite:///./molWH.db"
    SECRET_KEY:str = "bc0db4a2d1062ecaff76294bf8609fa33c4c401bac322e5bcbd07c5c5e048885"
    ALGORITHM :str = "HS256"
    ACCESS_TOKEN_EXPIRE_MINUTES: int = 30
    
    class Config:
        env_file = ".env"

@lru_cache()
def get_settings() -> Settings:
    return Settings()

settings  = get_settings()

engine = create_engine(settings.DB_URL, connect_args={"check_same_thread": False} \
                       if "sqlite" in settings.DB_URL else {},\
                        pool_pre_ping=True,\
                            echo = True)


