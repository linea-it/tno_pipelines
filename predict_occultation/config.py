from typing import Optional
from pydantic import BaseModel


class Config(BaseModel):
  task_id: int = 1
  job_id: Optional[int] = None
  asteroid_name: str = "2008 RH167"
  asteroid_path: str = "/app/predict_occultation/process001"


if __name__ == "__main__":
  import yaml

  cfg = Config()

  with open('config.yaml', 'w') as outfile:
    data_json = cfg.model_dump()
    print(data_json)
    yaml.dump(data_json, outfile)
