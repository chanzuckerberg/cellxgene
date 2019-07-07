import React from "react";
import { connect } from "react-redux";

import actions from "../../actions";

@connect(state => ({
	universe: state.universe,
	obsAnnotations: state.universe.obsAnnotations,
	saveInProgress: state.autosave?.saveInProgress ?? false,
	lastSavedObsAnnotations: state.autosave?.lastSavedObsAnnotations,
	error: state.autosave?.error,
	writableCategoriesEnabled: state.config?.parameters?.["label_file"] ?? false
}))
class Autosave extends React.Component {
	constructor(props) {
		super(props);
		this.state = {
			timer: null
		};
	}

	componentDidMount() {
		const { writableCategoriesEnabled } = this.props;

		let { timer } = this.state;
		if (timer) clearInterval(timer);
		if (writableCategoriesEnabled) {
			timer = setInterval(this.tick, 5000);
		} else {
			timer = null;
		}
		this.setState({ timer });
	}

	componentWillUnmount() {
		const { timer } = this.state;
		if (timer) this.clearInterval(timer);
	}

	tick = () => {
		const { dispatch, saveInProgress } = this.props;
		if (this.needToSave() && !saveInProgress) {
			dispatch(actions.saveObsAnnotations());
		}
	};

	needToSave = () => {
		/* return true if we need to save, false if we don't */
		const { obsAnnotations, lastSavedObsAnnotations } = this.props;
		return (
			lastSavedObsAnnotations && obsAnnotations !== lastSavedObsAnnotations
		);
	};

	statusMessage() {
		const { error } = this.props;
		if (error) {
			return `Autosave error: ${error}`;
		}
		return this.needToSave() ? "Unsaved" : "All saved";
	}

	render() {
		const { writableCategoriesEnabled } = this.props;
		return writableCategoriesEnabled ? (
			<div
				id="autosave"
				style={{
					position: "fixed",
					display: "inherit",
					right: 0,
					bottom: 0
				}}
			>
				{this.statusMessage()}
			</div>
		) : null;
	}
}

export default Autosave;
