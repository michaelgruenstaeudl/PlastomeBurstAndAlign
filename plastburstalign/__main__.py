from .user_parameters import UserParametersScript
from .plastome_burst_and_align import PlastomeRegionBurstAndAlign


if __name__ == "__main__":
    params = UserParametersScript()
    burst_align = PlastomeRegionBurstAndAlign(params)
    burst_align.execute()

    print("\nend of script\n")
