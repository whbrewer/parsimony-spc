%include('header')
<body>
%include('navbar')
%include('apps/alert')
<div class="container-fluid">
<form class="form-horizontal" action="/confirm" method="post" novalidate>
<input type="hidden" name="app" value="{{app}}">
%if defined('cid'):
<input type="hidden" name="cid" value="{{cid}}">
%end
<div class="col-xs-12" style="height:5px"></div>
<div class="col-xs-12 visible-xs" style="height:5px"></div>
<div class="form-group">
	<div class="col-xs-2">
		<button type="submit" class="btn btn-success"> <!-- pull-right -->
		Continue <em class="glyphicon glyphicon-forward"></em> </button>
	</div>
	<label for="desc" style="text-align:right" class="control-label col-sm-4 hidden-xs">
		<a href="#" data-toggle="tooltip" title="Separate labels by commas">Labels:</a></label>
	<div class="hidden-xs col-sm-6">
		<input type="text" id="desc" name="desc" class="form-control" style="width:100%"
			data-role="tagsinput" title="e.g. v2.5.1,bottleneck">
	</div>
</div>
<div class="col-xs-12" style="height:5px"></div>
<ul class="nav nav-pills" role="tablist">
	<li role="presentation" class="active">
		<a href="#basic" aria-controls="home" role="tab" data-toggle="tab">basic</a>
	</li>
</ul>
<div class="tab-content">
<div role="tabpanel" class="tab-pane fade in active" id="basic">
	<div class="form-group">
		<label for="species" class="control-label col-xs-6">
			species:</label>
		<div class="col-xs-12 col-sm-6">
			<input type="text" class="form-control" name="species" value="{{species}}"/>
		</div>
	</div>
</div>

</div>
</form>
%include('footer')
