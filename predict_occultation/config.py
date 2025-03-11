from pydantic import BaseModel


class Config(BaseModel):
  message: str = "Hello World"


if __name__ == "__main__":
  import yaml

  cfg = Config()

  with open('config.yaml', 'w') as outfile:
    data_json = cfg.model_dump()
    print(data_json)
    yaml.dump(data_json, outfile)
