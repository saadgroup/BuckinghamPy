from buckinghampy.buckinghampi import BuckinghamPi,sp
try:
    import ipywidgets as widgets
    from ipywidgets import HBox, VBox, Layout, Box
    from IPython.display import display, clear_output, Math, Markdown
except:
    pass



class BuckinghamPiGui(object):

    def __init__(self):
        import sys
        inJupyter = sys.argv[-1].endswith('json')
        if inJupyter == False:
            raise Exception("Cannot instantiate the class in a non jupyter cell!")

        self.continuousUpdate=True
        self.style = {'description_width': 'initial'}
        self.txt_box_layout = Layout(width='auto', height='32px')
        self.children_vbox = [] # variable list
        self.panels={}
        self.stylingTab()

        self.__display(self.panels)

        self.on_change()
        self.generateButtonPressed = False

    def __display(self,obj):
        for key in obj.keys():
            display(obj[key])

    def tabs_disiplay(self):
        self.tab = widgets.Tab(children=self.tabs)
        self.tab.set_title(0, 'results')
        self.tab.set_title(1, 'Numerical Setup')
        self.tab.set_title(2, 'style')

    def stylingTab(self):

        self.num_var = widgets.BoundedIntText(
            value=0,
            description='Number of Variables:',
            continuous_update=self.continuousUpdate,
            style=self.style
        )

        self.top_panel=HBox(children=[self.num_var])
        self.panels['top_panel']= self.top_panel

    def create_variable_Hbox(self,idx):

        setattr(self, 'var_name_{}'.format(idx),
                widgets.Textarea(
                    placeholder='var name',
                    description='Name:',
                    layout=self.txt_box_layout),
                )

        setattr(self, 'var_dimensions_{}'.format(idx),
                widgets.Textarea(
                    placeholder='var dimensions',
                    description='Dimensions:',
                    layout=self.txt_box_layout)
                )

        setattr(self, 'var_select_{}'.format(idx),
                widgets.Checkbox(
                    value=False,
                    description='Non-repeating')
                )

        box_layout = Layout(display='flex',
                            flex_flow='row',
                            align_items='flex-start',
                            # border='solid',
                            width='auto')
        items = [getattr(self, 'var_name_{}'.format(idx)), getattr(self, 'var_dimensions_{}'.format(idx)), getattr(self, 'var_select_{}'.format(idx))]
        box = Box(children=items, layout=box_layout)


        setattr(self, 'var_{}'.format(idx),box)


    def create_bottom_panel(self,change):
        if len(self.children_vbox):
            del self.children_vbox[-1]
        list_var_num = len(self.children_vbox)
        print("children v box list",self.children_vbox)
        print("list variable number ",list_var_num)
        print("variable number ",self.num_var.value)
        counter = list_var_num
        while counter < self.num_var.value:
            counter += 1
            self.create_variable_Hbox(counter)
            self.children_vbox.append(getattr(self, 'var_{}'.format(counter)))


        counter = list_var_num
        while counter > self.num_var.value:
            del self.children_vbox[-1]
            delattr(self,'var_{}'.format(counter))
            delattr(self, 'var_name_{}'.format(counter))
            delattr(self, 'var_dimensions_{}'.format(counter))
            delattr(self,'var_select_{}'.format(counter))
            counter -= 1

        self.generate_button_widget = widgets.Button(
                                                description='Generate',
                                                disabled=False,
                                                button_style='', # 'success', 'info', 'warning', 'danger' or ''
                                                tooltip='Generate',
                                                icon='check'
                                            )

        self.children_vbox.append(self.generate_button_widget)

        self.bottom_panel = \
            VBox(children=self.children_vbox)

        self.uncheck_chk_boxes()

        self.panels['bottom_panel'] = self.bottom_panel
        clear_output()
        self.__display(self.panels)
        self.output = widgets.Output()
        display(self.output)

        for idx in range(1,self.num_var.value+1):
            select_obj = getattr(self, 'var_select_{}'.format(idx))
            select_obj.observe(self.var_checkboxes_enable, names='value')

        self.generate_button_widget.on_click(self.generate_pressed)

    def uncheck_chk_boxes(self):
        all_chk_boxes = [getattr(self, 'var_select_{}'.format(idx)) for idx in range(1, self.num_var.value + 1)]
        chk_boxes_vals = [chk_box.value for chk_box in all_chk_boxes]
        for idx, chk_box in enumerate(all_chk_boxes):
            chk_box.value = False

    def change_visibility_select_box(self):
        all_chk_boxes = [getattr(self, 'var_select_{}'.format(idx)) for idx in range(1,self.num_var.value+1)]
        chk_boxes_vals = [chk_box.value for chk_box in all_chk_boxes]
        try:
            idx_true = chk_boxes_vals.index(True) + 1
            if idx_true:
                for idx, chk_box in enumerate(all_chk_boxes):
                    if idx + 1 != idx_true:
                        chk_box.disabled = not(chk_box.disabled)
        except:
            for idx, chk_box in enumerate(all_chk_boxes):
                chk_box.disabled = False

    def var_checkboxes_enable(self,change):
        self.change_visibility_select_box()

    def on_change(self):
        self.num_var.observe(self.create_bottom_panel,names='value')



    def collect_data(self):
        self.data={}
        var_num = self.num_var.value
        self.data['var_num'] = var_num
        self.data['vars'] ={}
        for idx in range(1,var_num+1):
            var_name = getattr(self, 'var_name_{}'.format(idx)).value
            var_dimensions = getattr(self, 'var_dimensions_{}'.format(idx)).value
            var_select = getattr(self, 'var_select_{}'.format(idx)).value

            self.data['vars'][var_name] = {'dimensions':var_dimensions,'non_repeating':var_select}

    def generate_solution(self):
        problem = BuckinghamPi()
        for varname in self.data['vars'].keys():
            problem.add_variable(name=varname, dimensions=self.data['vars'][varname]['dimensions'],
                                 non_repeating=self.data['vars'][varname]['non_repeating'])
        problem.generate_pi_terms()
        self.data['sol'] = problem.pi_terms

    def generate_pressed(self, *args):
        self.collect_data()

        try:
            self.generate_solution()
            with self.output:
                clear_output()
                self.print()
        except Exception as e:
            with self.output:
                clear_output()
                print(e)

    def print(self):
        for set_num, space in enumerate(self.data['sol']):
            latex_str= '\\text{Set }'
            latex_str+='{}: \\quad'.format(set_num+1)
            for num, term in enumerate(space):
                latex_str += '\\pi_{} = '.format(num+1)+sp.latex(term)
                latex_str += '\\quad'
            display(Math(latex_str))
            display(Markdown('---'))