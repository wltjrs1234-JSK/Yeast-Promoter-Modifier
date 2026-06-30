import sys
sys.path.append(".")
import data_loader

def test():
    for sys_id in ["YGR192C", "YAL003W"]:
        print(f"\nFetching {sys_id} promoter...")
        data = data_loader.get_promoter_data(sys_id)
        if data:
            print("Symbol:", data["symbol"])
            print("Systematic Name:", data["systematic_name"])
            print("Sequence Length:", len(data["seq"]))
            print("Sequence Preview:", data["seq"][:60] + "...")
            print("Sequence End:", data["seq"][-60:])
        else:
            print("Failed to fetch promoter data.")

if __name__ == "__main__":
    test()
