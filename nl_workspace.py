class Nl_workspace():
    def __init__(self, MY_ROOT, workspace_name, worksapce_detail=None):
        self.my_root = DA_ROOT
        self.workspace_name = workspace_name
        self.workspace_detail = worksapce_detail
        self.workspace_subdir = os.path.join(self.my_root, workspace_name)
        print('Workspace: ' + self.workspace_subdir)

        
 
