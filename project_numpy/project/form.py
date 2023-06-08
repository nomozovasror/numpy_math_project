from django import forms

class MyForm(forms.Form):
    a = forms.IntegerField(label='a')
    b = forms.IntegerField(label='b')